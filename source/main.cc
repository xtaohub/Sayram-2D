/*
 * File:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */


#include <iostream>
#include "Grid.h"
#include "D.h"
#include "BoundaryConditions.h"
#include "FVMSolver.h"

int main() {
    // Define grid parameters
    int nx = 10;
    int ny = 10;
    double dx = 1.0;
    double dy = 1.0;

    // Create grid object
    Grid grid(nx, ny, dx, dy);

    // Create diffusion coefficients object
    D diffusion(grid);

    // Create boundary conditions object
    BoundaryConditions boundaryConditions(grid, diffusion);

    // Create FVM solver object
    FVMSolver solver(grid, diffusion, boundaryConditions);

    // Define initial condition for f
    Eigen::MatrixXd f(nx, ny);
    f.setZero(); // Example: initialize with zeros

    // Define time step
    double dt = 0.1;

    // Time loop for solving
    for (int i = 0; i < 100; ++i) {
        // Solve using FVM solver
        solver.solve(f, dt);

        // Output or visualize f at each time step
        std::cout << "Time step: " << i << std::endl;
        std::cout << f << std::endl;
    }

    return 0;
}

