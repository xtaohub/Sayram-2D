/*
 * File:        FVMsolver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */


#include "FVMSolver.h"

FVMSolver::FVMSolver(const Grid& grid, const D& diffusion, const BoundaryConditions& boundaryConditions)
    : grid(grid), boundaryConditions(boundaryConditions) {
    int nx = grid.nx();
    int ny = grid.ny();

    // Resize and initialize the coefficient matrix M
    M.resize(nx * ny, nx * ny);
    M.setZero();

    // Resize and initialize the right-hand side vector R
    R.resize(nx * ny);
    R.setZero();
}

void FVMSolver::solve(Eigen::MatrixXd& f, const D& diffusion, double dt) {
    assembleSystem(f, diffusion, dt);

    // Solve the linear system Mf = R for f
    // f = M.fullPivLu().solve(R);

    // Apply boundary conditions
    // boundaryConditions.applyBoundaryConditions(f);
}

void FVMSolver::assembleSystem(Eigen::MatrixXd& f, const D& diffusion, double dt) {
    int nx = grid.nx();
    int ny = grid.ny();

    // Fill the coefficient matrix M and the right-hand side vector R
    // Example: for a simple 1D case
    for (int i = 0; i < nx * ny; ++i) {
        M(i, i) = 2.0 + dt * diffusion.getDap(dt, i % nx, i / nx); // Example: diagonal entries
        if (i > 0) {
            M(i, i - 1) = -dt * diffusion.getDap(dt, (i - 1) % nx, (i - 1) / nx); // Example: off-diagonal entries
        }
        if (i < nx * ny - 1) {
            M(i, i + 1) = -dt * diffusion.getDap(dt, (i + 1) % nx, (i + 1) / nx); // Example: off-diagonal entries
        }
        R(i) = f(i) * dt; // Example: R = f * dt
    }
}
