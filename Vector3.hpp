#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>
#include <utility>
#include <iostream>

// vector: position, also color (r,g,b) (extended from smallpt)
struct Vector3 {
  double x, y, z;
  Vector3() : x(0), y(0), z(0) {}
  Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

  // Copy & Move constructor
  Vector3(const Vector3& other) : x(other.x), y(other.y), z(other.z) {}
  Vector3(Vector3&& other) : x(std::exchange(other.x, 0)), y(std::exchange(other.y, 0)), z(std::exchange(other.z, 0)) {}

  // Copy & Move assignment
  Vector3& operator=(const Vector3& other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }
  Vector3& operator=(Vector3&& other) {
    x = std::move(other.x);
    y = std::move(other.y);
    z = std::move(other.z);
    return *this;
  }

  double magnitude() const { return std::sqrt(magnitudeSquared()); }
  double magnitudeSquared() const { return x * x + y * y + z * z; }

  inline Vector3 normalize() { return (*this) * (1.0 / magnitude()); }
  inline Vector3 normalized() const { return Vector3((*this) * (1.0 / magnitude())); }

  inline Vector3 reflect(const Vector3& normal) const { return *this - normal * 2.0 * dot(normal, *this); }

  inline Vector3 operator-() const { return Vector3(-x, -y, -z); }
  inline Vector3 operator+(const Vector3& rhs) const { return Vector3(x + rhs.x, y + rhs.y, z + rhs.z); }
  inline Vector3 operator-(const Vector3& rhs) const { return Vector3(x - rhs.x, y - rhs.y, z - rhs.z); }
  inline Vector3 operator+(double rhs) const { return Vector3(x + rhs, y + rhs, z + rhs); }
  inline Vector3 operator-(double rhs) const { return Vector3(x - rhs, y - rhs, z - rhs); }
  inline Vector3 operator*(double scalar) const { return Vector3(x * scalar, y * scalar, z * scalar); }

  friend double dot(const Vector3& lhs, const Vector3& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }
  friend Vector3 cross(const Vector3& lhs, const Vector3& rhs) {
    return Vector3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x);
  }
	
	friend std::ostream& operator<<(std::ostream& out, const Vector3& vector) {
    return out << vector.x << ' ' << vector.y << ' ' << vector.z;
  }
 
  inline double max() const { return std::fmax(std::fmax(x, y), z); }
  inline Vector3 onb(const Vector3& n) const {
    Vector3 u, w, v = n;
    if (n.z < -0.9999999) {
      u = Vector3(0.0, -1.0, 0.0);
      w = Vector3(-1.0, 0.0, 0.0);
    } else {
      const double a = 1.0 / (1.0 + n.z);
      const double b = -n.x * n.y * a;
      u = Vector3(1.0 - n.x * n.x * a, b, -n.x);
      w = Vector3(b, 1.0 - n.y * n.y * a, -n.y);
    }
    return Vector3(dot(*this, Vector3(u.x, v.x, w.x)), dot((*this), Vector3(u.y, v.y, w.y)),
                   dot(*this, Vector3(u.z, v.z, w.z)));
  }
};
#endif