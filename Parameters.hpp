#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <cmath>

const double PI = 2 * std::acos(0.0);

// parameters
constexpr int minPathLength = 3;  // avoid sampling direct illumination
constexpr int maxPathLength = 20;
constexpr double glossiness = 25.0;
constexpr int pixelWidth = 640;
constexpr int pixelHeight = 480;
constexpr int N_Init = 10000;
constexpr double largeStepProb = 0.3;

// scene independent constants
constexpr int numRNGsPerEvent = 2;
constexpr int maxEvents = maxPathLength + 1;
constexpr int numStatesSubpath = (maxEvents + 2) * numRNGsPerEvent;
constexpr int numStates = numStatesSubpath * 2;

#endif