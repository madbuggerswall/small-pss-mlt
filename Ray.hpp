#ifndef RAY_HPP
#define RAY_HPP

#include "Vector3.hpp"

// ray-sphere intersection (extended from smallpt)
struct Ray {
  Vector3 origin, direction;
  Ray(){};
  Ray(const Vector3& origin, const Vector3& direction) : origin(origin), direction(direction) {}
  // Copy & Move constructor
  Ray(const Ray& other) : origin(other.origin), direction(other.direction) {}
  Ray(Ray&& other) : origin(std::move(other.origin)), direction(std::move(other.direction)) {}

  // Copy assignment
  Ray& operator=(const Ray& other) {
    origin = other.origin;
    direction = other.direction;
    return *this;
  }

  // Move assignment
  Ray& operator=(Ray&& other) {
    direction = std::move(other.direction);
    origin = std::move(other.origin);
    return *this;
  }
};
#endif