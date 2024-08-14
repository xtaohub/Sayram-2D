#ifndef COMMON_H_
#define COMMON_H_

#define EIGEN_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "Eigen/Dense"
#include "Eigen/Sparse" 
#include "Eigen/Core"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

// constants
const double gPI = 3.141592653589793238462;
const double gD2R = gPI / 180.0; // convert degree to radian
const double gC = 1;
const double gE0 = 0.511875; // MeV
const double gME = gE0 / (gC * gC); 
const double gRE = 6371000;

using namespace std;  //bad practice, used here for convenience

#endif
