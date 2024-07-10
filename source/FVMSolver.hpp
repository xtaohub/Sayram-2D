/*
 * File:        FVMsolver.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
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

struct Alpha_K{ // The AlphaK matrix for each cell with four faces: east,west,north,south. 
  NTPFA_nodes e;
  NTPFA_nodes w;
  NTPFA_nodes n;
  NTPFA_nodes s;
};

// struct Unit_vector{
//   Eigen::Vector2d e;
//   Eigen::Vector2d w;
//   Eigen::Vector2d n;
//   Eigen::Vector2d s;
// };

class FVMSolver {
  public:
    FVMSolver(const Grid& g_in, const D& d_in, const BoundaryConditions& bc_in);

    void solve(double dt);
    const Eigen::MatrixXd& f() const { return f_; }

    void initial();

  private:
    const Grid& g;
    const D& d; 
    const BoundaryConditions& bc;

    SpMat M_;
    std::vector<T> M_coeffs_;

    Eigen::MatrixXd f_; 
    Eigen::VectorXd S_;

    Eigen::MatrixXd Id_; 

    Eigen::Matrix<Alpha_K, Eigen::Dynamic, Eigen::Dynamic> alpha_K_;

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

    double hdx_; 
    double hdy_; 

    void assemble_M();
    void construct_alpha_K();

};

#endif /* FVM_SOLVER_H */
