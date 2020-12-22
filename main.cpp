// smallpssmlt, primary sample space MLT by Toshiya Hachisuka
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
const int maxPathLength = 24;
const double glossiness = 25.0;
const int pixelWidth = 640;
const int pixelHeight = 480;
const int N_Init = 10000;
const double largeStepProb = 0.3;

// scene independent constants
const int numRNGsPerEvent = 2;
const int maxEvents = maxPathLength + 1;
const int numStatesSubpath = (maxEvents + 2) * numRNGsPerEvent;
const int numStates = numStatesSubpath * 2;

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
Vec VecCosine(const Vec& n, const double g, const double rnd1, const double rnd2) {
  const double temp1 = 2.0 * PI * rnd1, temp2 = pow(rnd2, 1.0 / (g + 1.0));
  const double s = std::sin(temp1), c = std::cos(temp1), t = std::sqrt(1.0 - temp2 * temp2);
  return Vec(s * t, temp2, c * t).onb(n);
}

// Scene: radius, position, color, material
Sphere sph[] = {Sphere(6.0, Vec(10, 70, 51.6), Vec(100., 100., 100.), LGHT),
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

bool intersect(const Ray& ray, double& t, int& id) {  // ray-sphere intrsect.
  int n = sizeof(sph) / sizeof(Sphere);
  double distance, inf = 1e20;
  t = inf;
  for (int i = 0; i < n; i++) {
    distance = sph[i].intersect(ray);
    if (distance < t) {
      t = distance;
      id = i;
    }
  }
  return t < inf;
}

// compensated summation
struct KahanAdder {
  double sum, carry, y;
  KahanAdder(const double b = 0.0) {
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
struct Camera {
  Vec origin, u, v, w;
  double dist;
  Camera() = default;
  Camera(const Vec& origin, const Vec& lookAt, const double fov) : origin(origin) {
    dist = pixelHeight / (2.0 * tan((fov / 2.0) * (PI / 180.0)));
    w = (lookAt - origin).norm();
    u = (w % Vec(0.0, 1.0, 0.0)).norm();
    v = u % w;
  }
  void set(const Vec& origin, const Vec& lookAt, const double fov) {
    this->origin = origin;
    dist = pixelHeight / (2.0 * tan((fov / 2.0) * (PI / 180.0)));
    w = (lookAt - origin).norm();
    u = (w % Vec(0.0, 1.0, 0.0)).norm();
    v = u % w;
  }
};

Camera camera;
const int camera_id = -2;
// image
Vec image[pixelWidth * pixelHeight];

// tone mapping
int toInt(double x) { return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5); }

// path data
struct Vertex {
  Vec point, normal;
  int id;
  Vertex(){};
  Vertex(const Vec& point, const Vec& normal, int id) : point(point), normal(normal), id(id) {}
};

struct Path {
  Vertex vertices[maxEvents];
  int vertexCount;
  Path() { vertexCount = 0; }

  Vertex& operator[](int index) { return vertices[index]; }
  Vertex operator[](int index) const { return vertices[index]; }
};

struct Contribution {
  double x, y;
  Vec color;
  Contribution(){};
  Contribution(double x, double y, const Vec& color) : x(x), y(y), color(color) {}
};

struct PathContribution {
 private:
  Contribution contributions[maxEvents * maxEvents];

 public:
  int contributionCount;
  double sc;
  PathContribution() : contributionCount(0), sc(0.0) {}

  Contribution& operator[](int index) { return contributions[index]; }
  Contribution operator[](int index) const { return contributions[index]; }
};

void AccumulatePathContribution(const PathContribution& pathContrib, const double mScaling) {
  if (pathContrib.sc == 0) return;
  for (int i = 0; i < pathContrib.contributionCount; i++) {
    const int ix = int(pathContrib[i].x), iy = int(pathContrib[i].y);
    const Vec color = pathContrib[i].color * mScaling;
    if ((ix < 0) || (ix >= pixelWidth) || (iy < 0) || (iy >= pixelHeight)) continue;
    image[ix + iy * pixelWidth] = image[ix + iy * pixelWidth] + color;
  }
}

// primary space Markov chain
inline double perturb(const double value, const double s1, const double s2) {
  double result;
  double randomValue = rnd();
  if (randomValue < 0.5) {
    randomValue = randomValue * 2.0;
    result = value + s2 * exp(-log(s2 / s1) * randomValue);
    if (result > 1.0) result -= 1.0;
  } else {
    randomValue = (randomValue - 0.5) * 2.0;
    result = value - s2 * exp(-log(s2 / s1) * randomValue);
    if (result < 0.0) result += 1.0;
  }
  return result;
}
struct MarkovChain {
  double states[numStates];
  PathContribution pathContribution;
  MarkovChain() {
    for (int i = 0; i < numStates; i++) states[i] = rnd();
  }
  MarkovChain largeStep() const {
    MarkovChain result;
    result.pathContribution = (*this).pathContribution;
    for (int i = 0; i < numStates; i++) result.states[i] = rnd();
    return result;
  }
  MarkovChain mutate() const {
    MarkovChain result;
    result.pathContribution = (*this).pathContribution;

    // pixel location
    result.states[0] = perturb(states[0], 2.0 / double(pixelWidth + pixelHeight), 0.1);
    result.states[1] = perturb(states[1], 2.0 / double(pixelWidth + pixelHeight), 0.1);

    // the rest
    for (int i = 2; i < numStates; i++) result.states[i] = perturb(states[i], 1.0 / 1024.0, 1.0 / 64.0);
    return result;
  }
};

// internal states of random numbers
int PathRndsOffset;
double prnds[numStates];
void InitRandomNumbersByChain(const MarkovChain& markovChain) {
  for (int i = 0; i < numStates; i++) prnds[i] = markovChain.states[i];
}
void InitRandomNumbers() {
  for (int i = 0; i < numStates; i++) prnds[i] = rnd();
}

// local sampling PDFs and standard terms
inline double GeometryTerm(const Vertex e0, const Vertex e1) {
  const Vec dv = e1.point - e0.point;
  const double d2 = dv.dot(dv);
  return fabs(e0.normal.dot(dv) * e1.normal.dot(dv)) / (d2 * d2);
}
inline double DirectionToArea(const Vertex current, const Vertex next) {
  const Vec dv = next.point - current.point;
  const double d2 = dv.dot(dv);
  return fabs(next.normal.dot(dv)) / (d2 * sqrt(d2));
}

inline double GlossyBRDF(const Vec& wi, const Vec& n, const Vec& wo) {
  const double won = wo.dot(n);
  const double win = wi.dot(n);
  const Vec r = (-wi).reflect(n);
  return (glossiness + 2.0) / (2.0 * PI) * pow(std::fmax(r.dot(wo), 0.0), glossiness) / std::fmax(fabs(win), fabs(won));
}
inline double GlossyPDF(const Vec& wi, const Vec& n, const Vec& wo) {
  const Vec r = (-wi).reflect(n);
  return (glossiness + 1.0) / (2.0 * PI) * pow(std::fmax(r.dot(wo), 0.0), glossiness);
}

inline double LambertianBRDF(const Vec& wi, const Vec& normal, const Vec& wo) { return 1.0 / PI; }
inline double LambertianPDF(const Vec& wi, const Vec& normal, const Vec& wo) { return fabs(wo.dot(normal)) / PI; }

// measurement contribution function
Vec PathThroughput(const Path& Xb) {
  Vec f = Vec(1.0, 1.0, 1.0);
  for (int i = 0; i < Xb.vertexCount; i++) {
    if (i == 0) {
      double W = 1.0 / double(pixelWidth * pixelHeight);
      Vec d0 = Xb[1].point - Xb[0].point;
      const double dist2 = d0.dot(d0);
      d0 = d0 * (1.0 / sqrt(dist2));
      const double c = d0.dot(camera.w);
      const double ds2 = (camera.dist / c) * (camera.dist / c);
      W = W / (c / ds2);
      f = f * (W * fabs(d0.dot(Xb[1].normal) / dist2));
    } else if (i == (Xb.vertexCount - 1)) {
      if (sph[Xb[i].id].refl == LGHT) {
        const Vec d0 = (Xb[i - 1].point - Xb[i].point).norm();
        const double L = LambertianBRDF(d0, Xb[i].normal, d0);
        f = f.mul(sph[Xb[i].id].color * L);
      } else {
        f = f * 0.0;
      }
    } else {
      const Vec d0 = (Xb[i - 1].point - Xb[i].point).norm();
      const Vec d1 = (Xb[i + 1].point - Xb[i].point).norm();
      double BRDF = 0.0;
      if (sph[Xb[i].id].refl == DIFF) {
        BRDF = LambertianBRDF(d0, Xb[i].normal, d1);
      } else if (sph[Xb[i].id].refl == GLOS) {
        BRDF = GlossyBRDF(d0, Xb[i].normal, d1);
      }
      f = f.mul(sph[Xb[i].id].color * BRDF * GeometryTerm(Xb[i], Xb[i + 1]));
    }
    if (f.Max() == 0.0) return f;
  }
  return f;
}

// check if the path can be connected or not (visibility term)
bool isConnectable(const Path& Xeye, const Path& Xlight, double& px, double& py) {
  Vec direction;
  const Vertex& Xeye_e = Xeye[Xeye.vertexCount - 1];
  const Vertex& Xlight_e = Xlight[Xlight.vertexCount - 1];

  bool result;
  if ((Xeye.vertexCount == 0) && (Xlight.vertexCount >= 2)) {
    // no direct hit to the film (pinhole)
    result = false;
  } else if ((Xeye.vertexCount >= 2) && (Xlight.vertexCount == 0)) {
    // direct hit to the light source
    result = (sph[Xeye_e.id].refl == LGHT);
    direction = (Xeye[1].point - Xeye[0].point).norm();
  } else if ((Xeye.vertexCount == 1) && (Xlight.vertexCount >= 1)) {
    // light tracing
    Ray ray(Xeye[0].point, (Xlight_e.point - Xeye[0].point).norm());
    double t;
    int id;
    result = (intersect(ray, t, id)) && (id == Xlight_e.id);
    direction = ray.direction;
  } else {
    // shadow ray connection
    Ray ray(Xeye_e.point, (Xlight_e.point - Xeye_e.point).norm());
    double t;
    int id;
    result = (intersect(ray, t, id)) && (id == Xlight_e.id);
    direction = (Xeye[1].point - Xeye[0].point).norm();
  }

  // get the pixel location
  Vec screenCenter = camera.origin + (camera.w * camera.dist);
  Vec screenPosition = camera.origin + (direction * (camera.dist / direction.dot(camera.w))) - screenCenter;
  px = camera.u.dot(screenPosition) + (pixelWidth * 0.5);
  py = -camera.v.dot(screenPosition) + (pixelHeight * 0.5);
  return result && ((px >= 0) && (px < pixelWidth) && (py >= 0) && (py < pixelHeight));
}

// path probability density
// - take the sum of all possible probability densities if the numbers of subpath vertices are not specified
double PathProbablityDensity(const Path& sampledPath, const int pathLength, const int specifiedNumEyeVertices = -1,
                             const int specifiedNumLightVertices = -1) {
  KahanAdder SumPDFs(0.0);
  bool specified = (specifiedNumEyeVertices != -1) && (specifiedNumLightVertices != -1);

  // number of eye subpath vertices
  for (int numEyeVertices = 0; numEyeVertices <= pathLength + 1; numEyeVertices++) {
    // extended BPT
    double p = 1.0;

    // number of light subpath vertices
    int numLightVertices = (pathLength + 1) - numEyeVertices;

    // we have pinhole camera
    if (numEyeVertices == 0) continue;

    // add all?
    if (specified && ((numEyeVertices != specifiedNumEyeVertices) || (numLightVertices != specifiedNumLightVertices)))
      continue;

    // sampling from the eye
    for (int i = -1; i <= numEyeVertices - 2; i++) {
      if (i == -1) {
        // PDF of sampling the camera position (the same delta function with the scaling 1.0 for all the PDFs - they
        // cancel out)
        p = p * 1.0;
      } else if (i == 0) {
        p = p * 1.0 / double(pixelWidth * pixelHeight);
        Vec Direction0 = (sampledPath[1].point - sampledPath[0].point).norm();
        double CosTheta = Direction0.dot(camera.w);
        double DistanceToScreen2 = camera.dist / CosTheta;
        DistanceToScreen2 = DistanceToScreen2 * DistanceToScreen2;
        p = p / (CosTheta / DistanceToScreen2);

        p = p * DirectionToArea(sampledPath[0], sampledPath[1]);
      } else {
        // PDF of sampling ith vertex
        Vec Direction0 = (sampledPath[i - 1].point - sampledPath[i].point).norm();
        Vec Direction1 = (sampledPath[i + 1].point - sampledPath[i].point).norm();

        if (sph[sampledPath[i].id].refl == DIFF) {
          p = p * LambertianPDF(Direction0, sampledPath[i].normal, Direction1);
        } else if (sph[sampledPath[i].id].refl == GLOS) {
          p = p * GlossyPDF(Direction0, sampledPath[i].normal, Direction1);
        }
        p = p * DirectionToArea(sampledPath[i], sampledPath[i + 1]);
      }
    }

    if (p != 0.0) {
      // sampling from the light source
      for (int i = -1; i <= numLightVertices - 2; i++) {
        if (i == -1) {
          // PDF of sampling the light position (assume area-based sampling)
          p = p * (1.0 / light_area);
        } else if (i == 0) {
          Vec Direction0 = (sampledPath[pathLength - 1].point - sampledPath[pathLength].point).norm();
          p = p * LambertianPDF(sampledPath[pathLength].normal, sampledPath[pathLength].normal, Direction0);
          p = p * DirectionToArea(sampledPath[pathLength], sampledPath[pathLength - 1]);
        } else {
          // PDF of sampling (PathLength - i)th vertex
          Vec Direction0 = (sampledPath[pathLength - (i - 1)].point - sampledPath[pathLength - i].point).norm();
          Vec Direction1 = (sampledPath[pathLength - (i + 1)].point - sampledPath[pathLength - i].point).norm();

          if (sph[sampledPath[pathLength - i].id].refl == DIFF) {
            p = p * LambertianPDF(Direction0, sampledPath[pathLength - i].normal, Direction1);
          } else if (sph[sampledPath[pathLength - i].id].refl == GLOS) {
            p = p * GlossyPDF(Direction0, sampledPath[pathLength - i].normal, Direction1);
          }
          p = p * DirectionToArea(sampledPath[pathLength - i], sampledPath[pathLength - (i + 1)]);
        }
      }
    }

    if (specified && (numEyeVertices == specifiedNumEyeVertices) && (numLightVertices == specifiedNumLightVertices))
      return p;

    // sum the probability density (use Kahan summation algorithm to reduce numerical issues)
    SumPDFs.add(p);
  }
  return SumPDFs.sum;
}

// path sampling
void TracePath(Path& path, const Ray& ray, const int rayLevel, const int maxRayLevel) {
  if (rayLevel >= maxRayLevel) return;
  double t;
  int id;
  if (!intersect(ray, t, id)) return;
  const Sphere& obj = sph[id];
  Vec intersectionPoint = ray.origin + ray.direction * t;
  Vec normal = (intersectionPoint - obj.position).norm();
  normal = normal.dot(ray.direction) < 0 ? normal : normal * -1;

  // set path data
  path[path.vertexCount] = Vertex(intersectionPoint, normal, id);
  path.vertexCount++;
  const double rnd0 = prnds[(rayLevel - 1) * numRNGsPerEvent + 0 + PathRndsOffset];
  const double rnd1 = prnds[(rayLevel - 1) * numRNGsPerEvent + 1 + PathRndsOffset];

  Ray scatteredRay;
  scatteredRay.origin = intersectionPoint;
  if (obj.refl == DIFF) {
    scatteredRay.direction = VecCosine(normal, 1.0, rnd0, rnd1);
    TracePath(path, scatteredRay, rayLevel + 1, maxRayLevel);
  } else if (obj.refl == GLOS) {
    scatteredRay.direction = VecCosine(ray.direction.reflect(normal), glossiness, rnd0, rnd1);
    TracePath(path, scatteredRay, rayLevel + 1, maxRayLevel);
  }
}

Ray SampleLightSources(const double rnd1, const double rnd2) {
  const Vec direction = VecRandom(rnd1, rnd2);
  return Ray(sph[light_id].position + (direction * sph[light_id].radius), direction);
}
Path GenerateLightPath(const int maxLightEvents) {
  Path lightPath;
  lightPath.vertexCount = 0;

  if (maxLightEvents == 0) return lightPath;
  for (int i = 0; i < maxEvents; i++) lightPath[i].id = -1;
  PathRndsOffset = numStatesSubpath;

  Ray ray = SampleLightSources(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  const Vec n = ray.direction;
  PathRndsOffset += numRNGsPerEvent;

  ray.direction = VecCosine(n, 1.0, prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  PathRndsOffset += numRNGsPerEvent;

  lightPath[0] = Vertex(ray.origin, n, light_id);
  lightPath.vertexCount++;
  TracePath(lightPath, ray, 1, maxLightEvents);
  return lightPath;
}

Ray SampleCamera(const double rnd1, const double rnd2) {
  const Vec su = camera.u * -(0.5 - rnd1) * pixelWidth;
  const Vec sv = camera.v * (0.5 - rnd2) * pixelHeight;
  const Vec sw = camera.w * camera.dist;
  return Ray(camera.origin, (su + sv + sw).norm());
}
Path GenerateEyePath(const int maxEyeEvents) {
  Path result;
  result.vertexCount = 0;

  if (maxEyeEvents == 0) return result;
  for (int i = 0; i < maxEvents; i++) result[i].id = -1;
  PathRndsOffset = 0;

  Ray ray = SampleCamera(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  PathRndsOffset += numRNGsPerEvent;

  result[0] = Vertex(ray.origin, camera.w, camera_id);
  result.vertexCount++;
  TracePath(result, ray, 1, maxEyeEvents);
  return result;
}

// balance heuristic
double MISWeight(const Path& sampledPath, const int numEyeVertices, const int numLightVertices, const int pathLength) {
  const double p_i = PathProbablityDensity(sampledPath, pathLength, numEyeVertices, numLightVertices);
  const double p_all = PathProbablityDensity(sampledPath, pathLength);
  if ((p_i == 0.0) || (p_all == 0.0)) {
    return 0.0;
  } else {
    return std::fmax(std::fmin(p_i / p_all, 1.0), 0.0);
  }
}

// BPT connections
// - limit the connection to a specific technique if s and t are provided
PathContribution CombinePaths(const Path& eyePath, const Path& lightPath, const int specifiedNumEyeVertices = -1,
                              const int specifiedNumLightVertices = -1) {
  PathContribution result;
  result.contributionCount = 0;
  result.sc = 0.0;
  const bool specified = (specifiedNumEyeVertices != -1) && (specifiedNumLightVertices != -1);

  // maxEvents = the maximum number of vertices
  for (int pathLength = minPathLength; pathLength <= maxPathLength; pathLength++) {
    for (int numEyeVertices = 0; numEyeVertices <= pathLength + 1; numEyeVertices++) {
      const int numLightVertices = (pathLength + 1) - numEyeVertices;

      if (numEyeVertices == 0) continue;  // no direct hit to the film (pinhole)
      if (numEyeVertices > eyePath.vertexCount) continue;
      if (numLightVertices > lightPath.vertexCount) continue;

      // take only the specified technique if provided
      if (specified && ((specifiedNumEyeVertices != numEyeVertices) || (specifiedNumLightVertices != numLightVertices)))
        continue;

      // extract subpaths
      Path Eyesubpath = eyePath;
      Path Lightsubpath = lightPath;
      Eyesubpath.vertexCount = numEyeVertices;
      Lightsubpath.vertexCount = numLightVertices;

      // check the path visibility
      double px = -1.0, py = -1.0;
      if (!isConnectable(Eyesubpath, Lightsubpath, px, py)) continue;

      // construct a full path
      Path sampledPath;
      for (int i = 0; i < numEyeVertices; i++) sampledPath[i] = eyePath[i];
      for (int i = 0; i < numLightVertices; i++) sampledPath[pathLength - i] = lightPath[i];
      sampledPath.vertexCount = numEyeVertices + numLightVertices;

      // evaluate the path
      Vec f = PathThroughput(sampledPath);
      double p = PathProbablityDensity(sampledPath, pathLength, numEyeVertices, numLightVertices);
      double w = MISWeight(sampledPath, numEyeVertices, numLightVertices, pathLength);
      if ((w <= 0.0) || (p <= 0.0)) continue;

      Vec c = f * (w / p);
      if (c.Max() <= 0.0) continue;

      // store the pixel contribution
      result[result.contributionCount] = Contribution(px, py, c);
      result.contributionCount++;

      // scalar contribution function
      result.sc = std::fmax(c.Max(), result.sc);

      // return immediately if the technique is specified
      if (specified && (specifiedNumEyeVertices == numEyeVertices) && (specifiedNumLightVertices == numLightVertices))
        return result;
    }
  }
  return result;
}

// main rendering process
int main(int argc, char* argv[]) {
  unsigned long samples = 0;
  const int ltime = (argc >= 2) ? std::fmax(atoi(argv[1]), 0) : 60 * 3;
  std::cout << ltime << std::endl;

  camera.set(Vec(50.0, 40.8, 220.0), Vec(50.0, 40.8, 0.0), 40.0);

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
  MarkovChain current, proposal;
  InitRandomNumbersByChain(current);
  current.pathContribution = CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents));

  // integration
  for (;;) {
    samples++;
    std::cout << "\rSamples: " << samples << std::flush;
    if (samples % pixelWidth) {
      gettimeofday(&currentTime, NULL);
      if (ltime < ((currentTime.tv_sec - startTime.tv_sec) + (currentTime.tv_usec - startTime.tv_usec) * 1.0E-6)) break;
    }

    // sample the path
    double isLargeStepDone;
    if (rnd() <= largeStepProb) {
      proposal = current.largeStep();
      isLargeStepDone = 1.0;
    } else {
      proposal = current.mutate();
      isLargeStepDone = 0.0;
    }
    InitRandomNumbersByChain(proposal);
    proposal.pathContribution = CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents));

    double a = 1.0;
    if (current.pathContribution.sc > 0.0)
      a = std::fmax(std::fmin(1.0, proposal.pathContribution.sc / current.pathContribution.sc), 0.0);

    // accumulate samples
    if (proposal.pathContribution.sc > 0.0)
      AccumulatePathContribution(proposal.pathContribution,
                                 (a + isLargeStepDone) / (proposal.pathContribution.sc / b + largeStepProb));
    if (current.pathContribution.sc > 0.0)
      AccumulatePathContribution(current.pathContribution,
                                 (1.0 - a) / (current.pathContribution.sc / b + largeStepProb));

    // update the chain
    if (rnd() <= a) current = proposal;
  }

  std::cout << std::endl << "Writing file..." << std::endl;
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