#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "Vector3.hpp"

enum Refl_t { DIFF, GLOS, LGHT };  // material types, used in radiance()

struct Sphere {
  double radius;
  Vec position, color;
  Refl_t refl;
  Sphere(double radius, Vec position, Vec color, Refl_t re_) :
      radius(radius),
      position(position),
      color(color),
      refl(re_) {}
  inline double intersect(const Ray& ray) const {  // returns distance
    Vec op = position - ray.o;
    double t, b = op.dot(ray.d), det = b * b - op.dot(op) + radius * radius;
    if (det < 0)
      return 1e20;
    else
      det = sqrt(det);
    return (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : 1e20);
  }
};
#endif