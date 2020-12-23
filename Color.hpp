#ifndef COLOR_HPP
#define COLOR_HPP

#include <cmath>

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
#endif