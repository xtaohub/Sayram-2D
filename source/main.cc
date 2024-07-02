/*
 * fin:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#define EIGEN_USE_THREADS
#include <iostream>
#include <emmintrin.h>
#include <omp.h>
#include <cassert>
#include "Grid.h"
#include "D.h"
#include "BoundaryConditions.h"
#include "Parameters.hpp"
#include "FVMSolver.hpp"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double>T;


int main() {

    Eigen::initParallel();

    #pragma omp parallel
    {
        Eigen::setNbThreads(Eigen::nbThreads());
    }

    // Create grid object
    Grid grid(nx, ny, dx, dy);

    // Create diffusion coefficients object
    D diffusion(grid);

    diffusion.constructD(0.0);

    // // Create boundary conditions object
    // BoundaryConditions boundaryConditions(grid, diffusion);

    // Create FVM solver object
    // FVMSolver solver(grid, diffusion);

    // Define initial condition for f
    Eigen::VectorXd f(nx * ny);

    // TODO BoundaryConditions modularization
    BoundaryConditions boundary(0, 0.0);

    FVMSolver fvmSolver(grid, diffusion, boundary);
    
    fvmSolver.initial(f);

    Eigen::MatrixXd M_Integration(nx*ny, nx*ny);
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(nx * ny, nx * ny);
    Eigen::MatrixXd M_inv(nx*ny, nx*ny);

    std::vector<T> M_coefficients;
    SpMat M(nx*ny, nx*ny);
    Eigen::VectorXd S_(nx*ny);

    double alpha_K[4][2][ny][nx];

    fvmSolver.constructAlpha_K(alpha_K);

    Eigen::MatrixXd L, U;
    string path;
    Eigen::SparseLU<SpMat> solver;

    // Time loop for solving
    for (int k = 0; k < steps; ++k) {
        // Solve using FVM solver
        // solver.solve(f, dt);
        
        S_ = Eigen::VectorXd::Zero(nx*ny);
        M_coefficients.clear();

        fvmSolver.timeForward(f, alpha_K, S_, M_coefficients);

        M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());
        
        M = M * (dt / (dx * dy)) + Id;
        S_ = S_ * (dt / (dx * dy)) + f;


        // Eigen::SimplicialLLT<SpMat> solver;
        // Eigen::VectorXd x = solver.compute(M).solve(S_);

        // Solving:
        // Eigen::SimplicialCholesky<SpMat> chol(M);  // performs a Cholesky factorization of A
        // Eigen::VectorXd x = chol.solve(S_);         // use the factorization to solve for the given right hand side

        solver.analyzePattern(M);
        solver.factorize(M);
        Eigen::VectorXd x = solver.solve(S_);

        f = x;

        int numThreads = Eigen::nbThreads();

        std::cout << "Using " << numThreads << " threads for Eigen" << std::endl;

        // Output or visualize f at each time step
        std::cout << "Time step: " << k << std::endl;
        
        if(k % printstep == printstep - 1){
            path = "../output/SMPPFV/smppfv" + std::to_string(int((k + 1) / printstep));

            ofstream outFile(path);
            assert(outFile);
    
            for (double value : f) {
                outFile << value << std::endl;
            }
            outFile.close();
        }
    }

    return 0;
}

