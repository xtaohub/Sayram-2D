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

    int ind2to1(int i, int j) const { 
      // map 2d indices to 1, column major
      return j*g.nx()+i; 
    }

    double vertex_f(int imin, int jmin) const { 
      // 
      // obtain flux at vertex by four-point inteprolation.
      // here (i,j) is the smallest indices among the four cells.
      // 
      //          | 
      //   (i,j)  |  (i,j+1)  
      //  --------v-----------
      //  (i+1,j) |  (i+1,j+1)
      //          |
      //
      return (f_(imin,jmin) + f_(imin,jmin+1) + f_(imin+1,jmin) + f_(imin+1,jmin+1)) / 4.0;
    }

    void assemble();
    void construct_alpha_K();

    void alpha_K(const Eigen::Matrix2d& Lambda_K, const Eigen::Vector2d& K, const Eigen::Vector2d& A, const Eigen::Vector2d& B, double* alpha_KA, double* alpha_KB); 

    // add coefficients to the Matrix M_ for inner grids (north)
    void coeff_M_add_n(int i, int j, double a, double p, double u1, double u2);

    // add coefficients to the Matrix M_ for inner grids (south)
    void coeff_M_add_s(int i, int j, double a, double p, double u3, double u4);

    // add coefficients to the Matrix M_ for inner grids (east)
    void coeff_M_add_e(int i, int j, double a, double p, double u4, double u1);

    // add coefficients to the Matrix M_ for inner grids (west)
    void coeff_M_add_w(int i, int j, double a, double p, double u2, double u3);

    // add coefficients, but the edge condition (east edge flux is 0, ignored)
    void coeff_add_n_edge(int i, int j, double a, double p, double u1, double u2);
    void coeff_add_s_edge(int i, int j, double a, double p, double u3, double u4);
    void coeff_add_w_edge(int i, int j, double a, double p, double u2, double u3);

    double calMuK(double a_K, double a_L) {
      if (a_K != 0 || a_L != 0){
        return abs(a_L) / (abs(a_K) + abs(a_L));
      } else {
        return 0.5;
      }
    }

    // double bsigma_inner_cells() { 
    //
    // }
    //
    // double bsigma_boundary_cells(){
    //
    // }

    double bsigma_plus(double bsigma){
        return (std::abs(bsigma) + std::abs(bsigma))/2.0;
    }

    double bsigma_minus(double bsigma){
        return (std::abs(bsigma) - std::abs(bsigma))/2.0;
    }

    double init_f(double a, double p){
      return exp(-(p2e(p, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (p * p);
    }

};

#endif /* FVM_SOLVER_H */
