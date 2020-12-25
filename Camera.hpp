#ifndef CAMERA_HPP
#define CAMERA_HPP

#include "Vector3.hpp"
#include "Parameters.hpp"

// pinhole camera
struct Camera {
  Vector3 origin, u, v, w;
  double dist;
  Camera() = default;
  Camera(const Vector3& origin, const Vector3& lookAt, const double fov) : origin(origin) {
    dist = pixelHeight / (2.0 * std::tan((fov / 2.0) * (PI / 180.0)));
    w = (lookAt - origin).normalize();
    u = cross(w, Vector3(0.0, 1.0, 0.0)).normalize();
    v = cross(u, w);
  }
};

#endif