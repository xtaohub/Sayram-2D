/*
 * fin:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
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
// #include <suitesparse/umfpack.h>


typedef Eigen::Triplet<double> T;

int main() {

    // Create grid object
    Grid g(nx, ny, dx, dy, dt);

    // Create diffusion coefficients object
    D diffusion(g);
    diffusion.constructD(0.0);

    // TODO BoundaryConditions modularization
    BoundaryConditions boundary(0, 0.0);

    FVMSolver solver(g, diffusion, boundary);
    
    solver.initial();

    string path;

    // The timer
    clock_t start, end;
    double cpu_time;
    start = clock();

    // Time loop for solving
    for (int k = 0; k < steps; ++k) {
        
        // Solve using FVM solver
        solver.update();

        // Output or visualize f at each time step
        std::cout << "Time step: " << k << std::endl;
        
        if(k % printstep == printstep - 1){
            path = "./output/SMPPFV/smppfv" + std::to_string(int((k + 1) / printstep));
    
            ofstream outFile(path);
            assert(outFile);
    
            for (double value : solver.f().reshaped()) {
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

