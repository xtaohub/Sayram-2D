/*
 * File:        BCs.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef BCS_H
#define BCS_H

#include "common.h"
#include "utils.h"
#include "Parameters.h"

class BCs {
public:
    BCs(const Parameters& p_in): paras(p_in) {};

    // Define your boundary condition functions here
    double init_f(double a0, double p) const{
      return exp(-(p2e(p, gE0) - 0.2) / 0.1) * (sin(a0) - sin(paras.alpha0_lc())) / (p * p);
    }

    double alpha0_lc(double t, double p) const{
      return 0.0;
    }

    double pmin(double t, double a0) const{
      return init_f(a0, paras.pmin()); 
    }

    double pmax(double t, double a0) const{
      return 0.0;
    }

private:
    const Parameters& paras; 

};

#endif /* BOUNDARY_CONDITIONS_H */
