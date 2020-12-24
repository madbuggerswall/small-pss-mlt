#ifndef COLOR_HPP
#define COLOR_HPP

#include <cmath>
#include <utility>

struct Color {
  double red, green, blue;
  Color() : red(0), green(0), blue(0) {}
  Color(double red, double green, double blue) : red(red), green(green), blue(blue) {}
 // Copy & Move constructor
  Color(const Color& other) : red(other.red), green(other.green), blue(other.blue) {}
  Color(Color&& other) : red(std::exchange(other.red, 0)), green(std::exchange(other.green, 0)), blue(std::exchange(other.blue, 0)) {}

  // Copy & Move assignment
  Color& operator=(const Color& other) {
    red = other.red;
    green = other.green;
    blue = other.blue;
    return *this;
  }
  Color& operator=(Color&& other) {
    red = std::move(other.red);
    green = std::move(other.green);
    blue = std::move(other.blue);
    return *this;
  }
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
  inline Color& operator+=(const Color& rhs) {
    red += rhs.red;
    green += rhs.green;
    blue += rhs.blue;
    return *this;
  }
  inline Color& operator-=(const Color& rhs) {
    red -= rhs.red;
    green -= rhs.green;
    blue -= rhs.blue;
    return *this;
  }
  
	inline double max() const { return std::fmax(std::fmax(red, green), blue); }
};
#endif