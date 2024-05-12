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

#include "Grid.h"
#include "D.h"
#include "BoundaryConditions.h"

class FVMSolver {
public:
    FVMSolver(const Grid& grid, const D& diffusion, const BoundaryConditions& boundaryConditions);

    // Solve the diffusion equation to update f
    void solve(Eigen::MatrixXd& f, const D& diffusion, double dt);

private:
    const Grid& grid;
    const BoundaryConditions& boundaryConditions;

    Eigen::MatrixXd M;
    Eigen::VectorXd R;

    // Assemble the coefficient matrix M and the right-hand side vector R
    void assembleSystem(Eigen::MatrixXd& f, const D& diffusion, double dt);
};

#endif /* FVM_SOLVER_H */
