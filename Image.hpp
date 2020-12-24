#ifndef IMAGE_HPP
#define IMAGE_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Color.hpp"

class Image {
 public:
  size_t height;
  size_t width;
  std::vector<Color> pixels;
  Image(size_t height, size_t width) : height(height), width(width) { pixels.reserve(height * width); }

  Color& operator[](int index) { return pixels[index]; }
  Color operator[](int index) const { return pixels[index]; }

  // tone mapping
  int toInt(double x) { return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5); }

  void writeToFile(std::string fileName, int samplesPerPixel);
};

void Image::writeToFile(std::string fileName, int samplesPerPixel) {
  std::cout << std::endl << "Writing file..." << std::endl;

  if (fileName.length() > 0)
    fileName += ".ppm";
  else
    fileName = "output.ppm";

  std::ofstream outputFile(fileName, std::ios::binary);

  outputFile << "P6" << std::endl;
  outputFile << width << "	" << height << std::endl;
  outputFile << "255" << std::endl;

  //	Divide the color total by the number of samples.
  double s = double(width * height) / double(samplesPerPixel);
  for (int i = 0; i < width * height; i++) {
    outputFile << static_cast<unsigned char>(toInt(pixels[i].red * s))
               << static_cast<unsigned char>(toInt(pixels[i].green * s))
               << static_cast<unsigned char>(toInt(pixels[i].blue * s));
  }
}

#endif