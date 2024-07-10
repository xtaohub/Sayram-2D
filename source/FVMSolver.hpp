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

typedef Eigen::Triplet<double> T;


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
    FVMSolver(const Grid& grid, const D& diffusion, const BoundaryConditions& boundaryConditions);

    // Solve the diffusion equation to update f
    void solve(Eigen::MatrixXd& f, const D& diffusion, double dt);

    void initial(Eigen::VectorXd& f);
    void timeForward(Eigen::VectorXd& f, Eigen::VectorXd &S_, std::vector<T> &M_coefficients);

    void construct_alpha_K();

  private:
    const Grid& grid;
    const BoundaryConditions& boundaryConditions;

    Eigen::MatrixXd M;
    Eigen::VectorXd R;

    Eigen::Matrix<Alpha_K, Eigen::Dynamic, Eigen::Dynamic> alpha_K_;

    Eigen::MatrixXd f_; 

    // Unit_vector nK_;
    // Unit_vector nT_; 

    double hdx_; 
    double hdy_; 
    // Assemble the coefficient matrix M and the right-hand side vector R
    void assembleSystem(Eigen::MatrixXd& f, const D& diffusion, double dt);
};

#endif /* FVM_SOLVER_H */
