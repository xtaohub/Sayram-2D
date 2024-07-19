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

#include "Mesh.h"
#include "D.h"
#include "BoundaryConditions.h"
#include "Eigen/Core"
#include "Parameters.hpp"
#include "Eigen/Sparse"
#include "Array.h"

struct NTPFA_nodes{ // two points A,B used in Nonlinear Two Point Approximation
  double A;
  double B;
}; 

class FVMSolver {
  public:
    FVMSolver(const Mesh& m_in, const D& d_in, const BoundaryConditions& bc_in);

    void update();
    const Eigen::MatrixXd& f() const { return f_; }

    void initial();

  private:
    const Mesh& m;
    const D& d; 
    const BoundaryConditions& bc;

    // M f = S
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

    SpMat V_;
    std::vector<T> V_coeffs_;
    SpMat M_;

    Eigen::MatrixXd f_;
    Eigen::MatrixXd U_; 
    Eigen::VectorXd S_;
    Eigen::VectorXd R_;

    Array<NTPFA_nodes, 3> alpha_K_; 

    //
    // use a matrix to store f at vertices to build a lookup 
    // table for fA and fB
    // vertex_f is of size (nx+1, ny+1)
    // 
    Array<double, 2> vertex_f_; 

    int ind2to1(int i, int j) const { // map 2d indices to 1, column major
      return j*m.nx()+i; 
    }

    void set_vertex_f(); 

    void assemble();
    void construct_alpha_K();

    void construct_U_();

    void alpha_K(const Eigen::Matrix2d& Lambda_K, const Eigen::Vector2d& K, const Eigen::Vector2d& A, const Eigen::Vector2d& B, double* alpha_KA, double* alpha_KB); 

    // add coefficient to M corresponds to the inbr cell of cell (i,j)
    // Here: the inbr neighbor is an inner cell.
    void coeff_M_add_inner(int i, int j, int inbr);  
  
    // Here: the inbr neighbor is a Dirichlet boundary cell.
    void coeff_M_add_dirbc(int i, int j, int inbr); 

    double coeff_a(double alphaA, double fA, double alphaB, double fB) const {
      return alphaA*fA + alphaB*fB; 
    }

    double coeff_mu(double aK, double aL) const {
      if (aK != 0 || aL != 0){
        return abs(aL) / (abs(aK) + abs(aL));
      } else {
        return 0.5;
      }
    }

    double bsigma_inner_cells(double aK, double muK, double aL, double muL) const { 
      return muL*aL-muK*aK;  
    }

    double bsigma_boundary_cells(double aK){
      return -aK; 
    }

    double bsigma_plus(double bsigma){
      return (std::abs(bsigma) + std::abs(bsigma))/2.0;
    }

    double bsigma_minus(double bsigma){
      return (std::abs(bsigma) - std::abs(bsigma))/2.0;
    }

    double init_f(double a, double p){
      return exp(-(p2e(p, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (p * p);
    }

    double boundary_imin(double p){
        return 0.0;
    }

    double boundary_jmin(double a){
      return exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
    }

    double boundary_jmax(double a){
      return 0.0;
    }

};

#endif /* FVM_SOLVER_H */
