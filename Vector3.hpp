#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>

// vector: position, also color (r,g,b) (extended from smallpt)
struct Vec {
  double x, y, z;
  Vec(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

  inline Vec operator-() const { return Vec(-x, -y, -z); }
  inline Vec operator+(const Vec& rhs) const { return Vec(x + rhs.x, y + rhs.y, z + rhs.z); }
  inline Vec operator-(const Vec& rhs) const { return Vec(x - rhs.x, y - rhs.y, z - rhs.z); }
  inline Vec operator+(double rhs) const { return Vec(x + rhs, y + rhs, z + rhs); }
  inline Vec operator-(double rhs) const { return Vec(x - rhs, y - rhs, z - rhs); }
  inline Vec operator*(double rhs) const { return Vec(x * rhs, y * rhs, z * rhs); }
  inline Vec mul(const Vec& rhs) const { return Vec(x * rhs.x, y * rhs.y, z * rhs.z); }
  inline Vec norm() { return (*this) * (1.0 / sqrt(x * x + y * y + z * z)); }
  inline double dot(const Vec& rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }
  Vec operator%(const Vec& rhs) const {
    return Vec(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x);
  }
  inline Vec reflect(const Vec& n) const { return (*this) - n * 2.0 * n.dot((*this)); }
  inline double Max() const { return std::fmax(std::fmax(x, y), z); }
  inline Vec onb(const Vec& n) const {
    Vec u, w, v = n;
    if (n.z < -0.9999999) {
      u = Vec(0.0, -1.0, 0.0);
      w = Vec(-1.0, 0.0, 0.0);
    } else {
      const double a = 1.0 / (1.0 + n.z);
      const double b = -n.x * n.y * a;
      u = Vec(1.0 - n.x * n.x * a, b, -n.x);
      w = Vec(b, 1.0 - n.y * n.y * a, -n.y);
    }
    return Vec((*this).dot(Vec(u.x, v.x, w.x)), (*this).dot(Vec(u.y, v.y, w.y)), (*this).dot(Vec(u.z, v.z, w.z)));
  }
};

// ray-sphere intersection (extended from smallpt)
struct Ray {
  Vec origin, direction;
  Ray(){};
  Ray(const Vec& origin, const Vec& direction) : origin(origin), direction(direction) {}
};

#endif