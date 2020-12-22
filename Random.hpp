#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

namespace Random {
  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.0, 1.0);

  inline float fraction() {
    return distribution(generator);
    // return rand() / (RAND_MAX + 1.0);
  }

  inline float range(float min, float max) {
    // return distribution(generator);
    return min + (max - min) * fraction();
  }

  inline float rangeInt(int min, int max) {
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
    // return min + (max - min) * fraction();
  }
}

#endif