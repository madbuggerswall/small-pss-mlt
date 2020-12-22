#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>
#include "Vector3.hpp"
#include "Math.hpp"

namespace Random {
  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.0, 1.0);

  inline float fraction() {
    return distribution(generator);
    // return rand() / (RAND_MAX + 1.0);
  }

  inline float range(float min, float max) {
    // return distribution(generator);
    return min + (max - min) * fraction();
  }

  inline float rangeInt(int min, int max) {
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
    // return min + (max - min) * fraction();
  }

  inline Vector3 vector() { return Vector3(fraction(), fraction(), fraction()); }

  inline Vector3 vectorRange(float min, float max) {
    return Vector3(range(min, max), range(min, max), range(min, max));
  }

  Vector3 unitVector() {
    auto angle = range(0, 2 * Math::pi);
    auto z = range(-1, 1);
    auto r = std::sqrt(1 - z * z);
    return Vector3(r * std::cos(angle), r * std::sin(angle), z);
  }

  Vector3 vectorInUnitSphere() {
    while (true) {
      auto vector = Random::vectorRange(-1, 1);
      if (vector.magnitudeSquared() >= 1) continue;
      return vector;
    }
  }

  Vector3 vectorInUnitDisk() {
    while (true) {
      auto vector = Vector3(range(-1, 1), range(-1, 1), 0);
      if (vector.magnitudeSquared() >= 1) continue;
      return vector;
    }
  }

  Vector3 vectorInHemiSphere(const Vector3& normal) {
    Vector3 inUnitSphere = vectorInUnitSphere();
    if (dot(inUnitSphere, normal) > 0.0)
      return inUnitSphere;
    else
      return -inUnitSphere;
  }

  inline Vector3 cosineDirection() {
    auto r1 = fraction();
    auto r2 = fraction();
    auto phi = 2 * Math::pi * r1;

    auto x = std::cos(phi) * std::sqrt(r2);
    auto y = std::sin(phi) * std::sqrt(r2);
    auto z = std::sqrt(1 - r2);

    return Vector3(x, y, z);
  }

}

#endif