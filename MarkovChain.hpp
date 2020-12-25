#ifndef MARKOVCHAIN_HPP
#define MARKOVCHAIN_HPP

#include <array>

#include "Color.hpp"
#include "Image.hpp"
#include "Parameters.hpp"
#include "Random.hpp"
#include "Vector3.hpp"

// path data
struct Vertex {
  Vector3 point, normal;
  int id;

  Vertex() = default;
  Vertex(const Vector3& point, const Vector3& normal, int id) : point(point), normal(normal), id(id) {}

  // Copy constructor
  Vertex(const Vertex& other) : point(other.point), normal(other.normal), id(other.id) {}

  // Move constructor
  Vertex(Vertex&& other) : point(std::move(other.point)), normal(std::move(other.normal)), id(std::move(other.id)) {}

  // Copy assignment
  Vertex& operator=(const Vertex& other) {
    point = other.point;
    normal = other.normal;
    id = other.id;
    return *this;
  }

  // Move assignment
  Vertex& operator=(Vertex&& other) {
    point = std::move(other.point);
    normal = std::move(other.normal);
    id = std::move(other.id);
    return *this;
  }
};

struct Path {
  std::array<Vertex, maxEvents> vertices;
  int vertexCount;

  Path() { vertexCount = 0; }

  // Copy constructor
  Path(const Path& other) : vertices(other.vertices), vertexCount(other.vertexCount) {}

  // Move constructor
  Path(Path&& other) : vertices(std::move(other.vertices)), vertexCount(std::move(other.vertexCount)) {}

  // Copy assignment
  Path& operator=(const Path& other) {
    vertices = other.vertices;
    vertexCount = other.vertexCount;
    return *this;
  }

  // Move assignment
  Path& operator=(Path&& other) {
    vertices = std::move(other.vertices);
    vertexCount = std::move(other.vertexCount);
    return *this;
  }

  Vertex& operator[](int index) { return vertices[index]; }
  Vertex operator[](int index) const { return vertices[index]; }
};

struct Contribution {
  double x, y;
  Color color;

  Contribution() {}

  Contribution(double x, double y, const Color& color) : x(x), y(y), color(color) {}

  // Copy constructor
  Contribution(const Contribution& other) : x(other.x), y(other.y), color(other.color) {}

  // Move constructor
  Contribution(Contribution&& other) : x(std::move(other.x)), y(std::move(other.y)), color(std::move(other.color)) {}

  // Copy assignment
  Contribution& operator=(const Contribution& other) {
    x = other.x;
    y = other.y;
    color = other.color;
    return *this;
  }

  // Move assignment
  Contribution& operator=(Contribution&& other) {
    x = std::move(other.x);
    y = std::move(other.y);
    color = std::move(other.color);
    return *this;
  }
};

struct PathContribution {
  std::array<Contribution, maxEvents * maxEvents> contributions;
  int contributionCount;
  double scalarContrib;

  PathContribution() : contributionCount(0), scalarContrib(0.0) {}

  // Copy constructor
  PathContribution(const PathContribution& other) :
      contributions(other.contributions),
      contributionCount(other.contributionCount),
      scalarContrib(other.scalarContrib) {}

  // Move constructor
  PathContribution(PathContribution&& other) :
      contributions(std::move(other.contributions)),
      contributionCount(std::move(other.contributionCount)),
      scalarContrib(std::move(other.scalarContrib)) {}

  // Copy assignment
  PathContribution& operator=(const PathContribution& other) {
    contributions = other.contributions;
    contributionCount = other.contributionCount;
    scalarContrib = other.scalarContrib;
    return *this;
  }

  // Move assignment
  PathContribution& operator=(PathContribution&& other) {
    contributions = std::move(other.contributions);
    // contributionCount = std::move(other.contributionCount);
    // scalarContrib = std::move(scalarContrib);
    contributionCount = other.contributionCount;
    scalarContrib = other.scalarContrib;
    return *this;
  }

  Contribution& operator[](int index) { return contributions[index]; }
  Contribution operator[](int index) const { return contributions[index]; }

  void accumulatePathContribution(const double scale, Image& image) {
    if (scalarContrib == 0) return;
    for (int i = 0; i < contributionCount; i++) {
      const int ix = int(contributions[i].x);
      const int iy = int(contributions[i].y);
      const Color color = contributions[i].color * scale;
      image[iy * pixelWidth + ix] += color;
    }
  }
};

struct MarkovChain {
  std::array<double, numStates> states;
  PathContribution pathContribution;

  MarkovChain() {
    for (int i = 0; i < numStates; i++) states[i] = Random::fraction();
  }

  // Copy constructor
  MarkovChain(const MarkovChain& other) : states(other.states), pathContribution(other.pathContribution) {}

  // Move constructor
  MarkovChain(MarkovChain&& other) :
      states(std::move(other.states)),
      pathContribution(std::move(other.pathContribution)) {}

  // Copy assignment
  MarkovChain& operator=(const MarkovChain& other) {
    states = other.states;
    pathContribution = other.pathContribution;
    return *this;
  }

  // Move assignment
  MarkovChain& operator=(MarkovChain&& other) {
    states = std::move(other.states);
    pathContribution = std::move(other.pathContribution);
    return *this;
  }

  // primary space Markov chain
  static inline double perturb(const double value, const double s1, const double s2) {
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

  MarkovChain largeStep() const {
    MarkovChain result;
    result.pathContribution = pathContribution;
    return result;
  }

  MarkovChain mutate() const {
    MarkovChain result;
    result.pathContribution = pathContribution;

    // pixel location
    result.states[0] = perturb(states[0], 2.0 / double(pixelWidth + pixelHeight), 0.1);
    result.states[1] = perturb(states[1], 2.0 / double(pixelWidth + pixelHeight), 0.1);

    // the rest
    for (int i = 2; i < numStates; i++) result.states[i] = perturb(states[i], 1.0 / 1024.0, 1.0 / 64.0);
    return result;
  }
};
#endif