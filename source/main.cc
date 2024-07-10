/*
 * fin:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#define EIGEN_USE_THREADS

// TODO MKL
// #define EIGEN_USE_MKL_ALL
// #define EIGEN_VECTORIZE_SSE4_2
#define EIGEN_NO_DEBUG
#include <iostream>
#include <emmintrin.h>
// TODO MKL
// #include "mkl.h"
#include <cassert>
#include "Grid.h"
#include "D.h"
#include "BoundaryConditions.h"
#include "Parameters.hpp"
#include "FVMSolver.hpp"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include <ctime>
// TODO: Have tried suitesparse, have little contribution of speed
// #include <suitesparse/umfpack.h>


typedef Eigen::Triplet<double> T;

int main() {

    // Create grid object
    Grid g(nx, ny, dx, dy);

    // Create diffusion coefficients object
    D diffusion(g);

    diffusion.constructD(0.0);

    // // Create boundary conditions object
    // BoundaryConditions boundaryConditions(g, diffusion);

    // Create FVM solver object
    // FVMSolver solver(g, diffusion);

    // Define initial condition for f

    // TODO BoundaryConditions modularization
    BoundaryConditions boundary(0, 0.0);

    FVMSolver solver(g, diffusion, boundary);
    
    solver.initial();
    solver.construct_alpha_K();

    string path;

    // The time test part
    clock_t start, end;
    double cpu_time;
    start = clock();

    // Time loop for solving
    for (int k = 0; k < steps; ++k) {
        // Solve using FVM solver
        // solver.solve(f, dt);

        solver.solve(dt);

        // Output or visualize f at each time step
        std::cout << "Time step: " << k << std::endl;
        
        if(k % printstep == printstep - 1){
            path = "./output/SMPPFV/smppfv" + std::to_string(int((k + 1) / printstep));
    
            ofstream outFile(path);
            assert(outFile);
    
            for (double value : solver.f()) {
                outFile << value << std::endl;
            }
            outFile.close();
        }
    }
    end = clock();
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "CPU time used " << cpu_time << " seconds" << std::endl;

    return 0;
}

