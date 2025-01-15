/*
 * File:        Albert_Young.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef ALBERT_YOUNG_H_
#define ALBERT_YOUNG_H_

#include "Equation.h"
#include "Albert_Young_IO.h"
#include "utils.h"


class Albert_Young:public Equation{
  public:
    Albert_Young(const Parameters& paras_in, const Mesh& m_in);

    double init_f(const Ind& ind) const; 
    void apply_bcs(Xtensor2d* vertex_fp) const; 
    void update(double t);  // where we update G, D, and Bcs.

    double xmin(double logE) const{
      return 0.0;
    }

    double ymin(double a0) const{
      return calculate_init_f(a0, paras.logEmin()); 
    }

    double ymax(double a0) const{
      return 0.0;
    }

  private:
    const Parameters& paras; 
    const Mesh& m;
    Albert_Young_IO io;

    // Define your boundary condition functions here
    double calculate_init_f(double a, double logE) const{
      double p = e2p(exp(logE), gE0);
      return exp(-(exp(logE) - 0.2) / 0.1) * (sin(a) - sin(5 * gPI / 180)) / (p * p) + gEPS;
    }

    double calculate_G(double alpha, double logE){
      double t = 1.30 - 0.56 * sin(alpha);
      return pow(e2p(exp(logE), gE0), 2) * t * sin(alpha) * cos(alpha) / dlogE_dp(logE, gE0);
    }

    void init(); 
    void locate(double alpha0, double logE, Loc* locp);
    void constructD(const Parameters& par, const Albert_Young_IO& io);
}; 


#endif /* ALBERT_YOUNG_H */

