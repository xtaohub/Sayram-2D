/*
 * File:        Equation.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef EQUATION_H_
#define EQUATION_H_

#include "common.h"
#include "utils.h"
#include "Mesh.h"
#include "Parameters.h"

//
// Definitions of the Equation are given here
// including D and boundary conditions.
// To handle different cases, we Equation as an abstract class.
//

class Equation{
  public:
    Equation(const Mesh& m) {
      std::size_t nx = m.nx();
      std::size_t ny = m.ny();

      G_.resize({nx,ny});

      Dxx_.resize({nx,ny});
      Dyy_.resize({nx,ny});
      Dxy_.resize({nx,ny});
    }

    double G(const Ind& ind) const { return G_(ind.i, ind.j); } 

    double Dxx(const Ind& ind) const { return Dxx_(ind.i, ind.j); }
    double Dyy(const Ind& ind) const { return Dyy_(ind.i, ind.j); } 
    double Dxy(const Ind& ind) const { return Dxy_(ind.i, ind.j); } 

    // pure virtual functions to be defined by the CASE of interest. 
    virtual double init_f(const Ind& ind) const = 0;
    virtual void apply_bcs(Xtensor2d* vertex_fp) const = 0; 
    virtual void update(double t) = 0;

  protected:
    Xtensor2d G_; 

    Xtensor2d Dxx_; 
    Xtensor2d Dyy_; 
    Xtensor2d Dxy_;
};


#endif /* EQUATION_H */

