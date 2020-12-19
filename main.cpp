#include <math.h>    // smallpssmlt, primary sample space MLT by Toshiya Hachisuka
#include <stdio.h>   // Usage: ./smallpssmlt time_sec
#include <stdlib.h>  // derived from smallpt, a path tracer by Kevin Beason, 2008

#include <cmath>  // derived from smallpt, a path tracer by Kevin Beason, 2008
#include <fstream>
#include <iostream>  // derived from smallpt, a path tracer by Kevin Beason, 2008

#include "Sphere.hpp"   // derived from smallpt, a path tracer by Kevin Beason, 2008
#include "Vector3.hpp"  // derived from smallpt, a path tracer by Kevin Beason, 2008
// 2015/01/26: Changed the default parameters. Fixed the bug that the normal was incorrectly flipped for diffuse
// surfaces (thanks to Hao Qin).

const double PI = 3.14159265358979;

// parameters
const int minPathLength = 3;  // avoid sampling direct illumination
const int maxPathLength = 20;
const double glossiness = 25.0;
const int pixelWidth = 640;
const int pixelHeight = 480;
const int N_Init = 10000;
const double largeStepProb = 0.3;

// scene independent constants
const int numRNGsPerEvent = 2;
const int maxEvents = (maxPathLength + 1);
const int numStatesSubpath = ((maxEvents + 2) * numRNGsPerEvent);
const int numStates = (numStatesSubpath * 2);

#include <sys/time.h>

// xorshift PRNG
double rnd(void) {
  static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 88675123;
  unsigned int t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8))) * (1.0 / 4294967296.0);
}

Vec VecRandom(const double rnd1, const double rnd2) {
  const double temp1 = 2.0 * PI * rnd1, temp2 = 2.0 * rnd2 - 1.0;
  const double s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2);
  return Vec(s * t, temp2, c * t);
}
Vec VecCosine(const Vec n, const double g, const double rnd1, const double rnd2) {
  const double temp1 = 2.0 * PI * rnd1, temp2 = pow(rnd2, 1.0 / (g + 1.0));
  const double s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2);
  return Vec(s * t, temp2, c * t).onb(n);
}

Sphere sph[] = {  // Scene: radius, position, color, material
    Sphere(6.0, Vec(10, 70, 51.6), Vec(100., 100., 100.), LGHT),
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(0.75, 0.25, 0.25), GLOS),
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(0.25, 0.25, 0.75), GLOS),
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(0.75, 0.65, 0.75), DIFF),
    Sphere(1e5, Vec(50, 40.8, -1e5 + 350), Vec(0.50, 0.50, 0.50), DIFF),
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(0.65, 0.75, 0.75), GLOS),
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(0.75, 0.75, 0.65), GLOS),
    Sphere(20, Vec(50, 20, 50), Vec(0.25, 0.75, 0.25), DIFF),
    Sphere(16.5, Vec(19, 16.5, 25), Vec(0.99, 0.99, 0.99), GLOS),
    Sphere(16.5, Vec(77, 16.5, 78), Vec(0.99, 0.99, 0.99), GLOS)};

const int light_id = 0;
const double light_area = (4.0 * PI * sph[light_id].radius * sph[light_id].radius);

bool intersect(const Ray& r, double& t, int& id) {  // ray-sphere intrsect.
  int n = sizeof(sph) / sizeof(Sphere);
  double d, inf = 1e20;
  t = inf;
  for (int i = 0; i < n; i++) {
    d = sph[i].intersect(r);
    if (d < t) {
      t = d;
      id = i;
    }
  }
  return t < inf;
}

// compensated summation
struct TKahanAdder {
  double sum, carry, y;
  TKahanAdder(const double b = 0.0) {
    sum = b;
    carry = 0.0;
    y = 0.0;
  }
  inline void add(const double b) {
    y = b - carry;
    const double t = sum + y;
    carry = (t - sum) - y;
    sum = t;
  }
};

// pinhole camera
struct TCamera {
  Vec o, u, v, w;
  double dist;
  void set(const Vec o_, const Vec l_, const double fov_) {
    o = o_;
    dist = pixelHeight / (2.0 * tan((fov_ / 2.0) * (PI / 180.0)));
    w = (l_ - o_).norm();
    u = (w % Vec(0.0, 1.0, 0.0)).norm();
    v = u % w;
  }
} Camera;

const int camera_id = -2;

// image
Vec image[pixelWidth * pixelHeight];

// tone mapping
int toInt(double x) { return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5); }

// path data
struct Vertex {
  Vec p, n;
  int id;
  Vertex(){};
  Vertex(Vec p_, Vec n_, int id_) : p(p_), n(n_), id(id_) {}
};
struct Path {
  Vertex vertices[maxEvents];
  int n;
  Path() { n = 0; }

  Vertex& operator[](int index) { return vertices[index]; }
  Vertex operator[](int index) const { return vertices[index]; }
};
struct Contribution {
  double x, y;
  Vec c;
  Contribution(){};
  Contribution(double x_, double y_, Vec c_) : x(x_), y(y_), c(c_) {}
};
struct PathContribution {
  Contribution c[maxEvents * maxEvents];
  int n;
  double sc;
  PathContribution() {
    n = 0;
    sc = 0.0;
  }
};
void AccumulatePathContribution(const PathContribution pc, const double mScaling) {
  if (pc.sc == 0) return;
  for (int i = 0; i < pc.n; i++) {
    const int ix = int(pc.c[i].x), iy = int(pc.c[i].y);
    const Vec c = pc.c[i].c * mScaling;
    if ((ix < 0) || (ix >= pixelWidth) || (iy < 0) || (iy >= pixelHeight)) continue;
    image[ix + iy * pixelWidth] = image[ix + iy * pixelWidth] + c;
  }
}

// primary space Markov chain
inline double perturb(const double value, const double s1, const double s2) {
  double Result;
  double r = rnd();
  if (r < 0.5) {
    r = r * 2.0;
    Result = value + s2 * exp(-log(s2 / s1) * r);
    if (Result > 1.0) Result -= 1.0;
  } else {
    r = (r - 0.5) * 2.0;
    Result = value - s2 * exp(-log(s2 / s1) * r);
    if (Result < 0.0) Result += 1.0;
  }
  return Result;
}
struct TMarkovChain {
  double u[numStates];
  PathContribution C;
  TMarkovChain() {
    for (int i = 0; i < numStates; i++) u[i] = rnd();
  }
  TMarkovChain large_step() const {
    TMarkovChain Result;
    Result.C = (*this).C;
    for (int i = 0; i < numStates; i++) Result.u[i] = rnd();
    return Result;
  }
  TMarkovChain mutate() const {
    TMarkovChain Result;
    Result.C = (*this).C;

    // pixel location
    Result.u[0] = perturb(u[0], 2.0 / double(pixelWidth + pixelHeight), 0.1);
    Result.u[1] = perturb(u[1], 2.0 / double(pixelWidth + pixelHeight), 0.1);

    // the rest
    for (int i = 2; i < numStates; i++) Result.u[i] = perturb(u[i], 1.0 / 1024.0, 1.0 / 64.0);
    return Result;
  }
};

// internal states of random numbers
int PathRndsOffset;
double prnds[numStates];
void InitRandomNumbersByChain(const TMarkovChain MC) {
  for (int i = 0; i < numStates; i++) prnds[i] = MC.u[i];
}
void InitRandomNumbers() {
  for (int i = 0; i < numStates; i++) prnds[i] = rnd();
}

// local sampling PDFs and standard terms
inline double GeometryTerm(const Vertex e0, const Vertex e1) {
  const Vec dv = e1.p - e0.p;
  const double d2 = dv.dot(dv);
  return fabs(e0.n.dot(dv) * e1.n.dot(dv)) / (d2 * d2);
}
inline double DirectionToArea(const Vertex current, const Vertex next) {
  const Vec dv = next.p - current.p;
  const double d2 = dv.dot(dv);
  return fabs(next.n.dot(dv)) / (d2 * sqrt(d2));
}

inline double GlossyBRDF(const Vec wi, const Vec n, const Vec wo) {
  const double won = wo.dot(n);
  const double win = wi.dot(n);
  const Vec r = (-wi).reflect(n);
  return (glossiness + 2.0) / (2.0 * PI) * pow(std::fmax(r.dot(wo), 0.0), glossiness) / std::fmax(fabs(win), fabs(won));
}
inline double GlossyPDF(const Vec wi, const Vec n, const Vec wo) {
  const Vec r = (-wi).reflect(n);
  return (glossiness + 1.0) / (2.0 * PI) * pow(std::fmax(r.dot(wo), 0.0), glossiness);
}

inline double LambertianBRDF(const Vec wi, const Vec n, const Vec wo) { return 1.0 / PI; }
inline double LambertianPDF(const Vec wi, const Vec n, const Vec wo) { return fabs(wo.dot(n)) / PI; }

// measurement contribution function
Vec PathThroughput(const Path Xb) {
  Vec f = Vec(1.0, 1.0, 1.0);
  for (int i = 0; i < Xb.n; i++) {
    if (i == 0) {
      double W = 1.0 / double(pixelWidth * pixelHeight);
      Vec d0 = Xb[1].p - Xb[0].p;
      const double dist2 = d0.dot(d0);
      d0 = d0 * (1.0 / sqrt(dist2));
      const double c = d0.dot(Camera.w);
      const double ds2 = (Camera.dist / c) * (Camera.dist / c);
      W = W / (c / ds2);
      f = f * (W * fabs(d0.dot(Xb[1].n) / dist2));
    } else if (i == (Xb.n - 1)) {
      if (sph[Xb[i].id].refl == LGHT) {
        const Vec d0 = (Xb[i - 1].p - Xb[i].p).norm();
        const double L = LambertianBRDF(d0, Xb[i].n, d0);
        f = f.mul(sph[Xb[i].id].color * L);
      } else {
        f = f * 0.0;
      }
    } else {
      const Vec d0 = (Xb[i - 1].p - Xb[i].p).norm();
      const Vec d1 = (Xb[i + 1].p - Xb[i].p).norm();
      double BRDF = 0.0;
      if (sph[Xb[i].id].refl == DIFF) {
        BRDF = LambertianBRDF(d0, Xb[i].n, d1);
      } else if (sph[Xb[i].id].refl == GLOS) {
        BRDF = GlossyBRDF(d0, Xb[i].n, d1);
      }
      f = f.mul(sph[Xb[i].id].color * BRDF * GeometryTerm(Xb[i], Xb[i + 1]));
    }
    if (f.Max() == 0.0) return f;
  }
  return f;
}

// check if the path can be connected or not (visibility term)
bool isConnectable(const Path Xeye, const Path Xlight, double& px, double& py) {
  Vec Direction;
  const Vertex& Xeye_e = Xeye[Xeye.n - 1];
  const Vertex& Xlight_e = Xlight[Xlight.n - 1];

  bool Result;
  if ((Xeye.n == 0) && (Xlight.n >= 2)) {
    // no direct hit to the film (pinhole)
    Result = false;
  } else if ((Xeye.n >= 2) && (Xlight.n == 0)) {
    // direct hit to the light source
    Result = (sph[Xeye_e.id].refl == LGHT);
    Direction = (Xeye[1].p - Xeye[0].p).norm();
  } else if ((Xeye.n == 1) && (Xlight.n >= 1)) {
    // light tracing
    Ray r(Xeye[0].p, (Xlight_e.p - Xeye[0].p).norm());
    double t;
    int id;
    Result = (intersect(r, t, id)) && (id == Xlight_e.id);
    Direction = r.d;
  } else {
    // shadow ray connection
    Ray r(Xeye_e.p, (Xlight_e.p - Xeye_e.p).norm());
    double t;
    int id;
    Result = (intersect(r, t, id)) && (id == Xlight_e.id);
    Direction = (Xeye[1].p - Xeye[0].p).norm();
  }

  // get the pixel location
  Vec ScreenCenter = Camera.o + (Camera.w * Camera.dist);
  Vec ScreenPosition = Camera.o + (Direction * (Camera.dist / Direction.dot(Camera.w))) - ScreenCenter;
  px = Camera.u.dot(ScreenPosition) + (pixelWidth * 0.5);
  py = -Camera.v.dot(ScreenPosition) + (pixelHeight * 0.5);
  return Result && ((px >= 0) && (px < pixelWidth) && (py >= 0) && (py < pixelHeight));
}

// path probability density
// - take the sum of all possible probability densities if the numbers of subpath vertices are not specified
double PathProbablityDensity(const Path SampledPath, const int PathLength, const int SpecifiedNumEyeVertices = -1,
                             const int SpecifiedNumLightVertices = -1) {
  TKahanAdder SumPDFs(0.0);
  bool Specified = (SpecifiedNumEyeVertices != -1) && (SpecifiedNumLightVertices != -1);

  // number of eye subpath vertices
  for (int NumEyeVertices = 0; NumEyeVertices <= PathLength + 1; NumEyeVertices++) {
    // extended BPT
    double p = 1.0;

    // number of light subpath vertices
    int NumLightVertices = (PathLength + 1) - NumEyeVertices;

    // we have pinhole camera
    if (NumEyeVertices == 0) continue;

    // add all?
    if (Specified && ((NumEyeVertices != SpecifiedNumEyeVertices) || (NumLightVertices != SpecifiedNumLightVertices)))
      continue;

    // sampling from the eye
    for (int i = -1; i <= NumEyeVertices - 2; i++) {
      if (i == -1) {
        // PDF of sampling the camera position (the same delta function with the scaling 1.0 for all the PDFs - they
        // cancel out)
        p = p * 1.0;
      } else if (i == 0) {
        p = p * 1.0 / double(pixelWidth * pixelHeight);
        Vec Direction0 = (SampledPath[1].p - SampledPath[0].p).norm();
        double CosTheta = Direction0.dot(Camera.w);
        double DistanceToScreen2 = Camera.dist / CosTheta;
        DistanceToScreen2 = DistanceToScreen2 * DistanceToScreen2;
        p = p / (CosTheta / DistanceToScreen2);

        p = p * DirectionToArea(SampledPath[0], SampledPath[1]);
      } else {
        // PDF of sampling ith vertex
        Vec Direction0 = (SampledPath[i - 1].p - SampledPath[i].p).norm();
        Vec Direction1 = (SampledPath[i + 1].p - SampledPath[i].p).norm();

        if (sph[SampledPath[i].id].refl == DIFF) {
          p = p * LambertianPDF(Direction0, SampledPath[i].n, Direction1);
        } else if (sph[SampledPath[i].id].refl == GLOS) {
          p = p * GlossyPDF(Direction0, SampledPath[i].n, Direction1);
        }
        p = p * DirectionToArea(SampledPath[i], SampledPath[i + 1]);
      }
    }

    if (p != 0.0) {
      // sampling from the light source
      for (int i = -1; i <= NumLightVertices - 2; i++) {
        if (i == -1) {
          // PDF of sampling the light position (assume area-based sampling)
          p = p * (1.0 / light_area);
        } else if (i == 0) {
          Vec Direction0 = (SampledPath[PathLength - 1].p - SampledPath[PathLength].p).norm();
          p = p * LambertianPDF(SampledPath[PathLength].n, SampledPath[PathLength].n, Direction0);
          p = p * DirectionToArea(SampledPath[PathLength], SampledPath[PathLength - 1]);
        } else {
          // PDF of sampling (PathLength - i)th vertex
          Vec Direction0 = (SampledPath[PathLength - (i - 1)].p - SampledPath[PathLength - i].p).norm();
          Vec Direction1 = (SampledPath[PathLength - (i + 1)].p - SampledPath[PathLength - i].p).norm();

          if (sph[SampledPath[PathLength - i].id].refl == DIFF) {
            p = p * LambertianPDF(Direction0, SampledPath[PathLength - i].n, Direction1);
          } else if (sph[SampledPath[PathLength - i].id].refl == GLOS) {
            p = p * GlossyPDF(Direction0, SampledPath[PathLength - i].n, Direction1);
          }
          p = p * DirectionToArea(SampledPath[PathLength - i], SampledPath[PathLength - (i + 1)]);
        }
      }
    }

    if (Specified && (NumEyeVertices == SpecifiedNumEyeVertices) && (NumLightVertices == SpecifiedNumLightVertices))
      return p;

    // sum the probability density (use Kahan summation algorithm to reduce numerical issues)
    SumPDFs.add(p);
  }
  return SumPDFs.sum;
}

// path sampling
void TracePath(Path& path, const Ray r, const int RayLevel, const int MaxRayLevel) {
  if (RayLevel >= MaxRayLevel) return;
  double t;
  int id;
  if (!intersect(r, t, id)) return;
  const Sphere& obj = sph[id];
  Vec p = r.o + r.d * t, n = (p - obj.position).norm();
  n = n.dot(r.d) < 0 ? n : n * -1;

  // set path data
  path[path.n] = Vertex(p, n, id);
  path.n++;
  const double rnd0 = prnds[(RayLevel - 1) * numRNGsPerEvent + 0 + PathRndsOffset];
  const double rnd1 = prnds[(RayLevel - 1) * numRNGsPerEvent + 1 + PathRndsOffset];

  Ray nr;
  nr.o = p;
  if (obj.refl == DIFF) {
    nr.d = VecCosine(n, 1.0, rnd0, rnd1);
    TracePath(path, nr, RayLevel + 1, MaxRayLevel);
  } else if (obj.refl == GLOS) {
    nr.d = VecCosine(r.d.reflect(n), glossiness, rnd0, rnd1);
    TracePath(path, nr, RayLevel + 1, MaxRayLevel);
  }
}

Ray SampleLightSources(const double rnd1, const double rnd2) {
  const Vec d = VecRandom(rnd1, rnd2);
  return Ray(sph[light_id].position + (d * sph[light_id].radius), d);
}
Path GenerateLightPath(const int MaxLightEvents) {
  Path Result;
  Result.n = 0;

  if (MaxLightEvents == 0) return Result;
  for (int i = 0; i < maxEvents; i++) Result[i].id = -1;
  PathRndsOffset = numStatesSubpath;

  Ray r = SampleLightSources(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  const Vec n = r.d;
  PathRndsOffset += numRNGsPerEvent;

  r.d = VecCosine(n, 1.0, prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  PathRndsOffset += numRNGsPerEvent;

  Result[0] = Vertex(r.o, n, light_id);
  Result.n++;
  TracePath(Result, r, 1, MaxLightEvents);
  return Result;
}

Ray SampleCamera(const double rnd1, const double rnd2) {
  const Vec su = Camera.u * -(0.5 - rnd1) * pixelWidth;
  const Vec sv = Camera.v * (0.5 - rnd2) * pixelHeight;
  const Vec sw = Camera.w * Camera.dist;
  return Ray(Camera.o, (su + sv + sw).norm());
}
Path GenerateEyePath(const int MaxEyeEvents) {
  Path Result;
  Result.n = 0;

  if (MaxEyeEvents == 0) return Result;
  for (int i = 0; i < maxEvents; i++) Result[i].id = -1;
  PathRndsOffset = 0;

  Ray r = SampleCamera(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  PathRndsOffset += numRNGsPerEvent;

  Result[0] = Vertex(r.o, Camera.w, camera_id);
  Result.n++;
  TracePath(Result, r, 1, MaxEyeEvents);
  return Result;
}

// balance heuristic
double MISWeight(const Path SampledPath, const int NumEyeVertices, const int NumLightVertices, const int PathLength) {
  const double p_i = PathProbablityDensity(SampledPath, PathLength, NumEyeVertices, NumLightVertices);
  const double p_all = PathProbablityDensity(SampledPath, PathLength);
  if ((p_i == 0.0) || (p_all == 0.0)) {
    return 0.0;
  } else {
    return std::fmax(std::fmin(p_i / p_all, 1.0), 0.0);
  }
}

// BPT connections
// - limit the connection to a specific technique if s and t are provided
PathContribution CombinePaths(const Path EyePath, const Path LightPath, const int SpecifiedNumEyeVertices = -1,
                              const int SpecifiedNumLightVertices = -1) {
  PathContribution Result;
  Result.n = 0;
  Result.sc = 0.0;
  const bool Specified = (SpecifiedNumEyeVertices != -1) && (SpecifiedNumLightVertices != -1);

  // maxEvents = the maximum number of vertices
  for (int PathLength = minPathLength; PathLength <= maxPathLength; PathLength++) {
    for (int NumEyeVertices = 0; NumEyeVertices <= PathLength + 1; NumEyeVertices++) {
      const int NumLightVertices = (PathLength + 1) - NumEyeVertices;

      if (NumEyeVertices == 0) continue;  // no direct hit to the film (pinhole)
      if (NumEyeVertices > EyePath.n) continue;
      if (NumLightVertices > LightPath.n) continue;

      // take only the specified technique if provided
      if (Specified && ((SpecifiedNumEyeVertices != NumEyeVertices) || (SpecifiedNumLightVertices != NumLightVertices)))
        continue;

      // extract subpaths
      Path Eyesubpath = EyePath;
      Path Lightsubpath = LightPath;
      Eyesubpath.n = NumEyeVertices;
      Lightsubpath.n = NumLightVertices;

      // check the path visibility
      double px = -1.0, py = -1.0;
      if (!isConnectable(Eyesubpath, Lightsubpath, px, py)) continue;

      // construct a full path
      Path SampledPath;
      for (int i = 0; i < NumEyeVertices; i++) SampledPath[i] = EyePath[i];
      for (int i = 0; i < NumLightVertices; i++) SampledPath[PathLength - i] = LightPath[i];
      SampledPath.n = NumEyeVertices + NumLightVertices;

      // evaluate the path
      Vec f = PathThroughput(SampledPath);
      double p = PathProbablityDensity(SampledPath, PathLength, NumEyeVertices, NumLightVertices);
      double w = MISWeight(SampledPath, NumEyeVertices, NumLightVertices, PathLength);
      if ((w <= 0.0) || (p <= 0.0)) continue;

      Vec c = f * (w / p);
      if (c.Max() <= 0.0) continue;

      // store the pixel contribution
      Result.c[Result.n] = Contribution(px, py, c);
      Result.n++;

      // scalar contribution function
      Result.sc = std::fmax(c.Max(), Result.sc);

      // return immediately if the technique is specified
      if (Specified && (SpecifiedNumEyeVertices == NumEyeVertices) && (SpecifiedNumLightVertices == NumLightVertices))
        return Result;
    }
  }
  return Result;
}

// main rendering process
int main(int argc, char* argv[]) {
  unsigned long samples = 0;
  const int ltime = (argc >= 2) ? std::fmax(atoi(argv[1]), 0) : 60 * 3;
  std::cout << ltime << std::endl;
  static const char* progr = "|*-/";
  FILE* f;

  Camera.set(Vec(50.0, 40.8, 220.0), Vec(50.0, 40.8, 0.0), 40.0);

  struct timeval startTime, currentTime;
  gettimeofday(&startTime, NULL);

  // PSSMLT
  // estimate normalization constant
  double b = 0.0;
  for (int i = 0; i < N_Init; i++) {
    fprintf(stderr, "\rPSSMLT Initializing: %5.2f", 100.0 * i / (N_Init));
    InitRandomNumbers();
    b += CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents)).sc;
  }
  b /= double(N_Init);
  fprintf(stderr, "\n");

  // initialize the Markov chain
  TMarkovChain current, proposal;
  InitRandomNumbersByChain(current);
  current.C = CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents));

  // integration
  for (;;) {
    samples++;
    fprintf(stderr, "\rPSSMLT %c", progr[(samples >> 8) & 3]);
    if (samples % pixelWidth) {
      gettimeofday(&currentTime, NULL);
      if (ltime < ((currentTime.tv_sec - startTime.tv_sec) + (currentTime.tv_usec - startTime.tv_usec) * 1.0E-6)) break;
    }

    // sample the path
    double isLargeStepDone;
    if (rnd() <= largeStepProb) {
      proposal = current.large_step();
      isLargeStepDone = 1.0;
    } else {
      proposal = current.mutate();
      isLargeStepDone = 0.0;
    }
    InitRandomNumbersByChain(proposal);
    proposal.C = CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents));

    double a = 1.0;
    if (current.C.sc > 0.0) a = std::fmax(std::fmin(1.0, proposal.C.sc / current.C.sc), 0.0);

    // accumulate samples
    if (proposal.C.sc > 0.0)
      AccumulatePathContribution(proposal.C, (a + isLargeStepDone) / (proposal.C.sc / b + largeStepProb));
    if (current.C.sc > 0.0) AccumulatePathContribution(current.C, (1.0 - a) / (current.C.sc / b + largeStepProb));

    // update the chain
    if (rnd() <= a) current = proposal;
  }

  // write out .ppm
  std::string fileName;
  if (argc > 2)
    fileName = std::string(argv[2]) + ".ppm";
  else
    fileName = "output.ppm";

  std::ofstream outputFile(fileName, std::ios::binary);

  outputFile << "P6" << std::endl;
  outputFile << pixelWidth << "	" << pixelHeight << std::endl;
  outputFile << "255" << std::endl;

 //	Divide the color total by the number of samples.
  double s = double(pixelWidth * pixelHeight) / double(samples);
  for (int i = 0; i < pixelWidth * pixelHeight; i++) {
    outputFile << static_cast<unsigned char>(toInt(image[i].x * s)) << static_cast<unsigned char>(toInt(image[i].y * s))
               << static_cast<unsigned char>(toInt(image[i].z * s));
  }
}