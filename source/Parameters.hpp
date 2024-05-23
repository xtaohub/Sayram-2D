#ifndef PARAMETERS_HPP_
#define PARAMETERS_HPP_

#include "utils.hpp"
#include "common.hpp"

const double ALPHA_LC = gPI * 5 / 180;
const double ALPHA_MAX = gPI / 2;
const double E_MIN = 0.2;
const double E_MAX = 5.0;
const double P_MIN = e2p(E_MIN, gE0);
const double P_MAX = e2p(E_MAX, gE0);
const int nx = 50;
const int ny = 50;
const double dx = (ALPHA_MAX - ALPHA_LC) / nx;
const double dy = (P_MAX - P_MIN) / ny;
const double hdx = dx / 2.0;
const double hdy = dy / 2.0;

const double E_RANGE[49] = {0.1000, 0.1085, 0.1177, 0.1277, 0.1385, 0.1503, 0.1631,
                           0.1769, 0.1919, 0.2082, 0.2259, 0.2451, 0.2659, 0.2885,
                           0.3130, 0.3396, 0.3684, 0.3997, 0.4336, 0.4704, 0.5104,
                           0.5537, 0.6008, 0.6518, 0.7071, 0.7671, 0.8323, 0.9030,
                           0.9796, 1.0630, 1.1530, 1.2510, 1.3570, 1.4720, 1.5970,
                           1.7330, 1.8800, 2.0400, 2.2130, 2.4010, 2.6050, 2.8260,
                           3.0660, 3.3270, 3.6090, 3.9150, 4.2480, 4.6090, 5.0000};

const double dlogE = (log(E_RANGE[48]) - log(E_RANGE[0])) / 48;

using namespace std; 


#endif /* PARAMETERS_H_ */