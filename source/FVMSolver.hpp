/*
 * File:        FVMsolver.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef FVM_SOLVER_H
#define FVM_SOLVER_H

// #include "mkl.h"
#include "Grid.h"
#include "D.h"
#include "BoundaryConditions.h"
#include "Eigen/Core"
#include "Parameters.hpp"
#include "Eigen/Sparse"


struct NTPFA_nodes{ // the two points A,B used in Nonlinear Two Point Approximation
  double A;
  double B;
}; 

struct Alpha_K{ // The AlphaK matrix for each cell with four faces: e(east) / w(west) corespond to alpha + da / - da faces; and n(north) / s(south) are p + dp / - dp faces. 
  NTPFA_nodes e;
  NTPFA_nodes w;
  NTPFA_nodes n;
  NTPFA_nodes s;
};

class FVMSolver {
  public:
    FVMSolver(const Grid& g_in, const D& d_in, const BoundaryConditions& bc_in);

    void update();
    const Eigen::MatrixXd& f() const { return f_; }

    void initial();

  private:
    const Grid& g;
    const D& d; 
    const BoundaryConditions& bc;

    // M f = S
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

    SpMat M_;
    std::vector<T> M_coeffs_;

    Eigen::MatrixXd f_; 
    Eigen::VectorXd S_;

    Eigen::MatrixXd Id_; 

    Eigen::Matrix<Alpha_K, Eigen::Dynamic, Eigen::Dynamic> alpha_K_;



    double hdx_; 
    double hdy_; 

    void assemble();
    void construct_alpha_K();

    double calMuK(double a_K, double a_L) {
      if (a_K != 0 || a_L != 0){
        return abs(a_L) / (abs(a_K) + abs(a_L));
      } else {
        return 0.5;
      }
    }

    double bsigma_inner_cells() { 

    }

    double bsigma_boundary_cells(){

    }

    double bsigma_plus(){
        return (std::abs(bsigma) + std::abs(bsigma))/2.0;
    }

    double bsigma_minus(){
        return (std::abs(bsigma) - std::abs(bsigma))/2.0;

    }

};

#endif /* FVM_SOLVER_H */
