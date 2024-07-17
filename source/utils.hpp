#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cmath>
#include <cassert>
#include "common.hpp"

inline double p2e(double p, double E0){ // convert momentum to energy
  return sqrt(p * p * gC * gC + E0 * E0) - E0; 
}

inline double e2p(double E, double E0){ // convert energy to momentum
  return sqrt(E * (E + 2 * E0)) / gC; 
}

inline double G(double alpha, double p){
  double t = 1.30 - 0.56 * sin(alpha);
  return p * p * t * sin(alpha) * cos(alpha);
}

#endif
