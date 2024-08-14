#ifndef UTILS_H_
#define UTILS_H_

#include <cmath>
#include "common.h"

inline double p2e(double p, double E0){ // convert momentum to energy
  return sqrt(p * p * gC * gC + E0 * E0) - E0; 
}

inline double e2p(double E, double E0){ // convert energy to momentum
  return sqrt(E * (E + 2 * E0)) / gC; 
}

#endif
