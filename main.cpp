// smallpssmlt, primary sample space MLT by Toshiya Hachisuka
// Usage: ./smallpssmlt time_sec
// derived from smallpt, a path tracer by Kevin Beason, 2008

// 2015/01/26: Changed the default parameters. Fixed the bug that the normal was incorrectly flipped for diffuse
// surfaces (thanks to Hao Qin).

#include <sys/time.h>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "Color.hpp"
#include "Image.hpp"
#include "Random.hpp"
#include "Ray.hpp"
#include "Sphere.hpp"
#include "Stopwatch.hpp"
#include "Vector3.hpp"

constexpr double PI = 3.14159265358979;

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
const int maxEvents = maxPathLength + 1;
const int numStatesSubpath = (maxEvents + 2) * numRNGsPerEvent;
const int numStates = numStatesSubpath * 2;

// Scene: radius, position, color, material
Sphere sph[] = {Sphere(6.0, Vector3(10, 70, 51.6), Color(100., 100., 100.), LGHT),
                Sphere(1e5, Vector3(1e5 + 1, 40.8, 81.6), Color(0.75, 0.25, 0.25), GLOS),
                Sphere(1e5, Vector3(-1e5 + 99, 40.8, 81.6), Color(0.25, 0.25, 0.75), GLOS),
                Sphere(1e5, Vector3(50, 40.8, 1e5), Color(0.75, 0.65, 0.75), DIFF),
                Sphere(1e5, Vector3(50, 40.8, -1e5 + 350), Color(0.50, 0.50, 0.50), DIFF),
                Sphere(1e5, Vector3(50, 1e5, 81.6), Color(0.65, 0.75, 0.75), GLOS),
                Sphere(1e5, Vector3(50, -1e5 + 81.6, 81.6), Color(0.75, 0.75, 0.65), GLOS),
                Sphere(20, Vector3(50, 20, 50), Color(0.25, 0.75, 0.25), DIFF),
                Sphere(16.5, Vector3(19, 16.5, 25), Color(0.99, 0.99, 0.99), GLOS),
                Sphere(16.5, Vector3(77, 16.5, 78), Color(0.99, 0.99, 0.99), GLOS)};

Vector3 VecRandom(const double rnd1, const double rnd2) {
  const double temp1 = 2.0 * PI * rnd1, temp2 = 2.0 * rnd2 - 1.0;
  const double s = std::sin(temp1), c = std::cos(temp1), t = std::sqrt(1.0 - temp2 * temp2);
  return Vector3(s * t, temp2, c * t);
}
Vector3 VecCosine(const Vector3& n, const double g, const double rnd1, const double rnd2) {
  const double temp1 = 2.0 * PI * rnd1, temp2 = pow(rnd2, 1.0 / (g + 1.0));
  const double s = std::sin(temp1), c = std::cos(temp1), t = std::sqrt(1.0 - temp2 * temp2);
  return Vector3(s * t, temp2, c * t).onb(n);
}

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
  Vector3 origin, u, v, w;
  double dist;
  Camera() = default;
  Camera(const Vector3& origin, const Vector3& lookAt, const double fov) : origin(origin) {
    dist = pixelHeight / (2.0 * std::tan((fov / 2.0) * (PI / 180.0)));
    w = (lookAt - origin).normalize();
    u = cross(w, Vector3(0.0, 1.0, 0.0)).normalize();
    v = cross(u, w);
  }
};

Camera camera;
const int camera_id = -2;
// image
Image image(pixelHeight, pixelWidth);

// path data
struct Vertex {
  Vector3 point, normal;
  int id;
  Vertex(){};
  Vertex(const Vector3& point, const Vector3& normal, int id) : point(point), normal(normal), id(id) {}
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
  Color color;
  Contribution(){};
  Contribution(double x, double y, const Color& color) : x(x), y(y), color(color) {}
};

struct PathContribution {
 private:
  Contribution contributions[maxEvents * maxEvents];

 public:
  int contributionCount;
  double scalarContrib;
  PathContribution() : contributionCount(0), scalarContrib(0.0) {}

  Contribution& operator[](int index) { return contributions[index]; }
  Contribution operator[](int index) const { return contributions[index]; }
};

void accumulatePathContribution(const PathContribution& pathContrib, const double scale) {
  if (pathContrib.scalarContrib == 0) return;
  for (int i = 0; i < pathContrib.contributionCount; i++) {
    const int ix = int(pathContrib[i].x);
    const int iy = int(pathContrib[i].y);
    const Color color = pathContrib[i].color * scale;
    if ((ix < 0) || (ix >= pixelWidth) || (iy < 0) || (iy >= pixelHeight)) {
      std::cout << "If ix&iy out of bounds." << std::endl;
      continue;
    }
    image[iy * pixelWidth + ix] += color;
  }
}

// primary space Markov chain
inline double perturb(const double value, const double s1, const double s2) {
  double result;
  double randomValue = Random::fraction();
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
    for (int i = 0; i < numStates; i++) states[i] = Random::fraction();
  }
  MarkovChain largeStep() const {
    MarkovChain result;
    result.pathContribution = (*this).pathContribution;
    for (int i = 0; i < numStates; i++) result.states[i] = Random::fraction();
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
  for (int i = 0; i < numStates; i++) prnds[i] = Random::fraction();
}

// local sampling PDFs and standard terms
inline double GeometryTerm(const Vertex& e0, const Vertex& e1) {
  const Vector3 dv = e1.point - e0.point;
  const double d2 = dot(dv, dv);
  return std::abs(dot(e0.normal, dv) * dot(e1.normal, dv)) / (d2 * d2);
}
inline double DirectionToArea(const Vertex& current, const Vertex& next) {
  const Vector3 dv = next.point - current.point;
  const double d2 = dot(dv, dv);
  return std::abs(dot(next.normal, dv)) / (d2 * sqrt(d2));
}

inline double GlossyBRDF(const Vector3& wi, const Vector3& n, const Vector3& wo) {
  const double won = dot(wo, n);
  const double win = dot(wi, n);
  const Vector3 reflected = (-wi).reflect(n);
  return (glossiness + 2.0) / (2.0 * PI) * pow(std::fmax(dot(reflected, wo), 0.0), glossiness) /
         std::fmax(std::abs(win), std::abs(won));
}
inline double GlossyPDF(const Vector3& wi, const Vector3& n, const Vector3& wo) {
  const Vector3 reflected = (-wi).reflect(n);
  return (glossiness + 1.0) / (2.0 * PI) * pow(std::fmax(dot(reflected, wo), 0.0), glossiness);
}

inline double LambertianBRDF(const Vector3& wi, const Vector3& normal, const Vector3& wo) { return 1.0 / PI; }
inline double LambertianPDF(const Vector3& wi, const Vector3& normal, const Vector3& wo) {
  return std::abs(dot(wo, normal)) / PI;
}

// measurement contribution function
Color PathThroughput(const Path& Xb) {
  Color color = Color(1.0, 1.0, 1.0);
  for (int i = 0; i < Xb.vertexCount; i++) {
    if (i == 0) {
      double W = 1.0 / double(pixelWidth * pixelHeight);
      Vector3 d0 = Xb[1].point - Xb[0].point;
      const double dist2 = dot(d0, d0);
      d0 = d0 * (1.0 / sqrt(dist2));
      const double c = dot(d0, camera.w);
      const double ds2 = (camera.dist / c) * (camera.dist / c);
      W /= (c / ds2);
      color *= (W * std::abs(dot(d0, Xb[1].normal) / dist2));
    } else if (i == (Xb.vertexCount - 1)) {
      if (sph[Xb[i].id].refl == LGHT) {
        const Vector3 d0 = (Xb[i - 1].point - Xb[i].point).normalize();
        const double L = LambertianBRDF(d0, Xb[i].normal, d0);
        color *= sph[Xb[i].id].color * L;
      } else {
        color *= 0.0;
      }
    } else {
      const Vector3 d0 = (Xb[i - 1].point - Xb[i].point).normalize();
      const Vector3 d1 = (Xb[i + 1].point - Xb[i].point).normalize();
      double BRDF = 0.0;
      if (sph[Xb[i].id].refl == DIFF) {
        BRDF = LambertianBRDF(d0, Xb[i].normal, d1);
      } else if (sph[Xb[i].id].refl == GLOS) {
        BRDF = GlossyBRDF(d0, Xb[i].normal, d1);
      }
      color *= sph[Xb[i].id].color * BRDF * GeometryTerm(Xb[i], Xb[i + 1]);
    }
    if (color.max() == 0.0) return color;
  }
  return color;
}

// check if the path can be connected or not (visibility term)
bool isConnectable(const Path& Xeye, const Path& Xlight, double& px, double& py) {
  Vector3 direction;
  const Vertex& Xeye_e = Xeye[Xeye.vertexCount - 1];
  const Vertex& Xlight_e = Xlight[Xlight.vertexCount - 1];

  bool result;
  if ((Xeye.vertexCount == 0) && (Xlight.vertexCount >= 2)) {
    // no direct hit to the film (pinhole)
    result = false;
  } else if ((Xeye.vertexCount >= 2) && (Xlight.vertexCount == 0)) {
    // direct hit to the light source
    result = (sph[Xeye_e.id].refl == LGHT);
    direction = (Xeye[1].point - Xeye[0].point).normalize();
  } else if ((Xeye.vertexCount == 1) && (Xlight.vertexCount >= 1)) {
    // light tracing
    Ray ray(Xeye[0].point, (Xlight_e.point - Xeye[0].point).normalize());
    double t;
    int id;
    result = (intersect(ray, t, id)) && (id == Xlight_e.id);
    direction = ray.direction;
  } else {
    // shadow ray connection
    Ray ray(Xeye_e.point, (Xlight_e.point - Xeye_e.point).normalize());
    double t;
    int id;
    result = (intersect(ray, t, id)) && (id == Xlight_e.id);
    direction = (Xeye[1].point - Xeye[0].point).normalize();
  }

  // get the pixel location
  Vector3 screenCenter = camera.origin + (camera.w * camera.dist);
  Vector3 screenPosition = camera.origin + (direction * (camera.dist / dot(direction, camera.w))) - screenCenter;
  px = dot(camera.u, screenPosition) + (pixelWidth * 0.5);
  py = dot(-camera.v, screenPosition) + (pixelHeight * 0.5);
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
    double pdfValue = 1.0;

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
        pdfValue *= 1.0;
      } else if (i == 0) {
        pdfValue *= 1.0 / double(pixelWidth * pixelHeight);
        Vector3 Direction0 = (sampledPath[1].point - sampledPath[0].point).normalize();
        double CosTheta = dot(Direction0, camera.w);
        double DistanceToScreen2 = camera.dist / CosTheta;
        DistanceToScreen2 = DistanceToScreen2 * DistanceToScreen2;
        pdfValue /= (CosTheta / DistanceToScreen2);

        pdfValue *= DirectionToArea(sampledPath[0], sampledPath[1]);
      } else {
        // PDF of sampling ith vertex
        Vector3 Direction0 = (sampledPath[i - 1].point - sampledPath[i].point).normalize();
        Vector3 Direction1 = (sampledPath[i + 1].point - sampledPath[i].point).normalize();

        if (sph[sampledPath[i].id].refl == DIFF) {
          pdfValue *= LambertianPDF(Direction0, sampledPath[i].normal, Direction1);
        } else if (sph[sampledPath[i].id].refl == GLOS) {
          pdfValue *= GlossyPDF(Direction0, sampledPath[i].normal, Direction1);
        }
        pdfValue *= DirectionToArea(sampledPath[i], sampledPath[i + 1]);
      }
    }

    if (pdfValue != 0.0) {
      // sampling from the light source
      for (int i = -1; i <= numLightVertices - 2; i++) {
        if (i == -1) {
          // PDF of sampling the light position (assume area-based sampling)
          pdfValue *= (1.0 / light_area);
        } else if (i == 0) {
          Vector3 Direction0 = (sampledPath[pathLength - 1].point - sampledPath[pathLength].point).normalize();
          pdfValue *= LambertianPDF(sampledPath[pathLength].normal, sampledPath[pathLength].normal, Direction0);
          pdfValue *= DirectionToArea(sampledPath[pathLength], sampledPath[pathLength - 1]);
        } else {
          // PDF of sampling (PathLength - i)th vertex
          Vector3 dir0 = (sampledPath[pathLength - (i - 1)].point - sampledPath[pathLength - i].point).normalize();
          Vector3 dir1 = (sampledPath[pathLength - (i + 1)].point - sampledPath[pathLength - i].point).normalize();

          if (sph[sampledPath[pathLength - i].id].refl == DIFF) {
            pdfValue *= LambertianPDF(dir0, sampledPath[pathLength - i].normal, dir1);
          } else if (sph[sampledPath[pathLength - i].id].refl == GLOS) {
            pdfValue *= GlossyPDF(dir0, sampledPath[pathLength - i].normal, dir1);
          }
          pdfValue *= DirectionToArea(sampledPath[pathLength - i], sampledPath[pathLength - (i + 1)]);
        }
      }
    }

    if (specified && (numEyeVertices == specifiedNumEyeVertices) && (numLightVertices == specifiedNumLightVertices))
      return pdfValue;

    // sum the probability density (use Kahan summation algorithm to reduce numerical issues)
    SumPDFs.add(pdfValue);
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
  Vector3 intersectionPoint = ray.origin + ray.direction * t;
  Vector3 normal = (intersectionPoint - obj.position).normalize();
  normal = dot(normal, ray.direction) < 0 ? normal : normal * -1;

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

// Get an initial ray coming out from a light source to a random direction.
Ray SampleLightSources(const double rnd1, const double rnd2) {
  const Vector3 direction = VecRandom(rnd1, rnd2);
  return Ray(sph[light_id].position + (direction * sph[light_id].radius), direction);
}

// Generates a path originating from a light source.
// maxLightEvents: Maximum number of vertices in an eye path.
Path GenerateLightPath(const int maxLightEvents) {
  Path lightPath;
  lightPath.vertexCount = 0;

  if (maxLightEvents == 0) return lightPath;
  for (int i = 0; i < maxEvents; i++) lightPath[i].id = -1;
  PathRndsOffset = numStatesSubpath;

  Ray ray = SampleLightSources(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  const Vector3 n = ray.direction;
  PathRndsOffset += numRNGsPerEvent;

  ray.direction = VecCosine(n, 1.0, prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
  PathRndsOffset += numRNGsPerEvent;

  lightPath[0] = Vertex(ray.origin, n, light_id);
  lightPath.vertexCount++;
  TracePath(lightPath, ray, 1, maxLightEvents);
  return lightPath;
}

// Get an initial ray coming out from the camera to a random direction.
Ray SampleCamera(const double rnd1, const double rnd2) {
  const Vector3 su = camera.u * -(0.5 - rnd1) * pixelWidth;
  const Vector3 sv = camera.v * (0.5 - rnd2) * pixelHeight;
  const Vector3 sw = camera.w * camera.dist;
  return Ray(camera.origin, (su + sv + sw).normalize());
}

// Generates a path originating from the camera.
// maxEyeEvents: Maximum number of vertices in an eye path.
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
  result.scalarContrib = 0.0;
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
      Color color = PathThroughput(sampledPath);
      double pathPDF = PathProbablityDensity(sampledPath, pathLength, numEyeVertices, numLightVertices);
      double weight = MISWeight(sampledPath, numEyeVertices, numLightVertices, pathLength);
      if ((weight <= 0.0) || (pathPDF <= 0.0)) continue;

      color *= (weight / pathPDF);
      if (color.max() <= 0.0) continue;

      // store the pixel contribution
      result[result.contributionCount] = Contribution(px, py, color);
      result.contributionCount++;

      // scalar contribution function
      result.scalarContrib = std::fmax(color.max(), result.scalarContrib);

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
  const int runningTime = (argc >= 2) ? std::fmax(atoi(argv[1]), 0) : 60 * 3;
  std::cout << "Run for " << runningTime << " seconds." << std::endl;

  camera = Camera(Vector3(50.0, 40.8, 220.0), Vector3(50.0, 40.8, 0.0), 40.0);

  Stopwatch stopwatch;

  // PSSMLT
  // estimate normalization constant
  double b = 0.0;
  for (int i = 0; i < N_Init; i++) {
    std::cout << std::fixed;
    std::cout << std::setprecision(1) << "\rPSSMLT Initializing: " << 100.0 * i / (N_Init) << std::flush;
    InitRandomNumbers();
    b += CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents)).scalarContrib;
  }
  b /= double(N_Init);
  fprintf(stderr, "\n");

  // initialize the Markov chain
  MarkovChain current, proposal;
  InitRandomNumbersByChain(current);
  current.pathContribution = CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents));

  // integration
  stopwatch.start();

  for (;;) {
    samples++;
    std::cout << "\rSamples: " << samples << std::flush;
    if (samples % pixelWidth) {
      if (stopwatch.getTime() > runningTime) break;
    }

    // sample the path
    double isLargeStepDone;
    if (Random::fraction() <= largeStepProb) {
      proposal = current.largeStep();
      isLargeStepDone = 1.0;
    } else {
      proposal = current.mutate();
      isLargeStepDone = 0.0;
    }
    InitRandomNumbersByChain(proposal);
    proposal.pathContribution = CombinePaths(GenerateEyePath(maxEvents), GenerateLightPath(maxEvents));

    double a = 1.0;
    if (current.pathContribution.scalarContrib > 0.0)
      a = std::fmax(std::fmin(1.0, proposal.pathContribution.scalarContrib / current.pathContribution.scalarContrib),
                    0.0);

    // accumulate samples
    if (proposal.pathContribution.scalarContrib > 0.0)
      accumulatePathContribution(proposal.pathContribution,
                                 (a + isLargeStepDone) / (proposal.pathContribution.scalarContrib / b + largeStepProb));
    if (current.pathContribution.scalarContrib > 0.0)
      accumulatePathContribution(current.pathContribution,
                                 (1.0 - a) / (current.pathContribution.scalarContrib / b + largeStepProb));

    // update the chain
    if (Random::fraction() <= a) current = proposal;
  }
  image.writeToFile(argv[2], samples);
}