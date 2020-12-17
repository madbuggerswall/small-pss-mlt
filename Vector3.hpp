#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>

// vector: position, also color (r,g,b) (extended from smallpt)
struct Vec {
  double x, y, z;
  Vec(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

  inline Vec operator-() const { return Vec(-x, -y, -z); }
  inline Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
  inline Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
  inline Vec operator+(double b) const { return Vec(x + b, y + b, z + b); }
  inline Vec operator-(double b) const { return Vec(x - b, y - b, z - b); }
  inline Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
  inline Vec mul(const Vec& b) const { return Vec(x * b.x, y * b.y, z * b.z); }
  inline Vec norm() { return (*this) * (1.0 / sqrt(x * x + y * y + z * z)); }
  inline double dot(const Vec& b) const { return x * b.x + y * b.y + z * b.z; }
  Vec operator%(const Vec& b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
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
  Vec o, d;
  Ray(){};
  Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

#endif