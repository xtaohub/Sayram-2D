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
    double init_f(double a, double p) const{
      return exp(-(p2e(p, gE0) - 0.2) / 0.1) * sin(a) / (p * p);
    }

    double pmin(double a) const{
      return init_f(a, paras.pmin()); 
    }

    double pmax(double a) const{
      return 0.0;
    }

private:
    const Parameters& paras; 

};

#endif /* BOUNDARY_CONDITIONS_H */
