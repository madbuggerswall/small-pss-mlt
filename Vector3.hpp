#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>

// vector: position, also color (r,g,b) (extended from smallpt)
struct Vector3 {
  double x, y, z;
  Vector3() : x(0), y(0), z(0) {}
  Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

  inline Vector3 operator-() const { return Vector3(-x, -y, -z); }
  inline Vector3 operator+(const Vector3& rhs) const { return Vector3(x + rhs.x, y + rhs.y, z + rhs.z); }
  inline Vector3 operator-(const Vector3& rhs) const { return Vector3(x - rhs.x, y - rhs.y, z - rhs.z); }
  inline Vector3 operator+(double rhs) const { return Vector3(x + rhs, y + rhs, z + rhs); }
  inline Vector3 operator-(double rhs) const { return Vector3(x - rhs, y - rhs, z - rhs); }
  inline Vector3 operator*(double scalar) const { return Vector3(x * scalar, y * scalar, z * scalar); }
  inline Vector3 mul(const Vector3& rhs) const { return Vector3(x * rhs.x, y * rhs.y, z * rhs.z); }
  inline Vector3 norm() { return (*this) * (1.0 / sqrt(x * x + y * y + z * z)); }
  inline double dot(const Vector3& rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }
  Vector3 operator%(const Vector3& rhs) const {
    return Vector3(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x);
  }
  inline Vector3 reflect(const Vector3& n) const { return (*this) - n * 2.0 * n.dot((*this)); }
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
    return Vector3((*this).dot(Vector3(u.x, v.x, w.x)), (*this).dot(Vector3(u.y, v.y, w.y)), (*this).dot(Vector3(u.z, v.z, w.z)));
  }
};

struct Color {
  double red, green, blue;
  Color() : red(0), green(0), blue(0) {}
  Color(double red, double green, double blue) : red(red), green(green), blue(blue) {}

  inline Color operator*(double scalar) const { return Color(red * scalar, green * scalar, blue * scalar); }
  inline Color operator*(const Color& rhs) const { return Color(red * rhs.red, green * rhs.green, blue * rhs.blue); }
  inline Color operator+(const Color& rhs) const { return Color(red + rhs.red, green + rhs.green, blue + rhs.blue); }
  inline Color operator-(const Color& rhs) const { return Color(red - rhs.red, green - rhs.green, blue - rhs.blue); }

  inline Color& operator*=(double scalar) {
    red *= scalar;
    green *= scalar;
    blue *= scalar;
    return *this;
  }
  inline Color& operator*=(const Color& rhs) {
    red *= rhs.red;
    green *= rhs.green;
    blue *= rhs.blue;
    return *this;
  }

  inline double max() const { return std::fmax(std::fmax(red, green), blue); }
};

// ray-sphere intersection (extended from smallpt)
struct Ray {
  Vector3 origin, direction;
  Ray(){};
  Ray(const Vector3& origin, const Vector3& direction) : origin(origin), direction(direction) {}
};

#endif