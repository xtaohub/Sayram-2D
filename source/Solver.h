/*
 * File:        Solver.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "common.h"
#include "Mesh.h"
#include "D.h"
#include "BCs.h"
#include "Parameters.h"
#include <vector>
#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"


struct NTPFA_node{ // two points A,B used in Nonlinear Two Point Approximation
  double A;
  double B;
}; 

class Solver {
  public:
    Solver(const Parameters& paras_in, const Mesh& m_in, const D& d_in, const BCs& bcs_in);

    void update();
    double t() const { return t_; }
    const Eigen::MatrixXd& f() const { return f_; }

  private:
    const Parameters& paras; 
    const Mesh& m;
    const D& d; 
    const BCs& bcs;

    double t_; 

    // M f = R
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

    SpMat M_;
    std::vector<T> M_coeffs_;

    Eigen::MatrixXd f_;
    Eigen::VectorXd R_;

    Eigen::MatrixXd tau_; 

    xt::xtensor<NTPFA_node,3> alpha_osf_; // alpha_one_sided_flux

    //
    // use a matrix to store f at vertices to build a lookup 
    // table for fA and fB
    // vertex_f is of size (nx+1, ny+1)
    // 
    xt::xtensor<double,2> vertex_f_; 

    void update_vertex_f(); 

    void assemble();

    void construct_alpha_osf();
    void alpha_osf_func(const Eigen::Matrix2d& Lambda_K, const Point& K, const Point& A, const Point& B, NTPFA_node* nodep);

    // add coefficients to M and R corresponds to the inbr cell of cell (i,j)
    // Here: the inbr neighbor is an inner cell.
    void coeff_add_inner(int i, int j, int inbr);  
  
    // Here: the inbr neighbor is a Dirichlet boundary cell.
    void coeff_add_dirbc(int i, int j, int inbr); 

    double coeff_mu(double aK, double aL) const {
      if (aK != 0 || aL != 0){
        return abs(aL) / (abs(aK) + abs(aL));
      } else {
        return 0.5;
      }
    }

    double bsigma_plus(double bsigma){
      return (std::abs(bsigma) + bsigma)/2.0;
    }

    double bsigma_minus(double bsigma){
      return (std::abs(bsigma) - bsigma)/2.0;
    }

    double G(double alpha, double p){ // this is the Jacobian for (a0, log(p))
      double T = 1.30 - 0.56 * sin(alpha);
      return p * p * p * T * sin(alpha) * cos(alpha);
    }

    void init();

    double bounce_period(double a, double p) const; 

};

#endif /* SOLVER_H */
