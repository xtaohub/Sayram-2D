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


typedef Eigen::SparseMatrix<double> SpMat;
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
    Eigen::VectorXd f(nx * ny);

    // TODO BoundaryConditions modularization
    BoundaryConditions boundary(0, 0.0);

    FVMSolver fvmSolver(g, diffusion, boundary);
    
    fvmSolver.initial(f);

    Eigen::MatrixXd M_Integration(nx*ny, nx*ny);
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(nx * ny, nx * ny);
    Eigen::MatrixXd M_inv(nx*ny, nx*ny);

    std::vector<T> M_coefficients;
    SpMat M(nx*ny, nx*ny);
    Eigen::VectorXd S_(nx*ny);


    // fvmSolver.construct_alpha_K();

    Eigen::MatrixXd L, U;
    string path;
    // Eigen::SparseLU<SpMat> solver;
    // TODO suitesparse part
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    Eigen::VectorXd x;

    // The time test part
    clock_t start, end;
    double cpu_time;
    start = clock();

    // Time loop for solving
    for (int k = 0; k < steps; ++k) {
        // Solve using FVM solver
        // solver.solve(f, dt);
        
        S_ = Eigen::VectorXd::Zero(nx*ny);
        M_coefficients.clear();

        fvmSolver.timeForward(f, S_, M_coefficients);

        M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());
        
        M = M * (dt / (dx * dy)) + Id;
        S_ = S_ * (dt / (dx * dy)) + f;

        // suitesparse part
        M.makeCompressed();
        
        solver.analyzePattern(M);
        solver.factorize(M);
        x = solver.solve(S_);
        //
        // f = x;

        // Output or visualize f at each time step
        // std::cout << "Time step: " << k << std::endl;
        
    //     if(k % printstep == printstep - 1){
    //         path = "./output/SMPPFV/smppfv" + std::to_string(int((k + 1) / printstep));
    //
    //         ofstream outFile(path);
    //         assert(outFile);
    // 
    //         for (double value : f) {
    //             outFile << value << std::endl;
    //         }
    //         outFile.close();
    //     }
    }
    end = clock();
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "CPU time used " << cpu_time << " seconds" << std::endl;

    return 0;
}

