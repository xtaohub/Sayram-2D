/*
 * fin:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */


#include <iostream>
#include <cassert>
#include "Grid.h"
#include "D.h"
#include "BoundaryConditions.h"
#include "Parameters.hpp"

int main() {
    // Create grid object
    Grid grid(nx, ny, dx, dy);

    // Create diffusion coefficients object
    D diffusion(grid);

    diffusion.constructD(0.0);

    // Create boundary conditions object
    // BoundaryConditions boundaryConditions(grid, diffusion);

    // // Create FVM solver object
    // FVMSolver solver(grid, diffusion, boundaryConditions);

    // Define initial condition for f
    Eigen::VectorXd f(nx * ny);
    for (int j = 0; j< ny; j++){
        for (int i = 0; i < nx; i++){
            double p = P_MIN + hdy + dy * j;
            double a = ALPHA_LC + hdx + dx * i;
            f(j * 50 + i) = exp(-(p2e(p, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (p * p);
        }
    }

    Eigen::MatrixXd M_(nx*ny, nx*ny);
    Eigen::VectorXd S_(nx*ny);

    // Define time step
    double dt = 10.0;

    double alpha_K[4][2][ny][nx];

    for (int j = 0; j< ny; j++){
        for (int i = 0; i < nx; i++){
            // 8 direction of vector:
            double p = P_MIN + hdy + dy * j;
            double a = ALPHA_LC + hdx + dx * i;

            Eigen::Matrix2d Lambda_K{{diffusion.getDaa(0.0, j, i), diffusion.getDap(0.0, j, i)}, {diffusion.getDap(0.0, j, i), diffusion.getDpp(0.0, j, i)}};
            Lambda_K *= G(a, p);
            
            // 0, 1 represent counterclockwise
            // 1 represents +1
            // e, w, n, s represent East, West, North, South edges
            
            // East:
            double sigma_e = dy;
            Eigen::MatrixXd n_Ke{{1, 0}}; 
            Eigen::MatrixXd n_1_Te{{hdx, hdy}};
            Eigen::MatrixXd n_1_Re{{-hdy}, {hdx}};
            Eigen::MatrixXd n_0_Te{{hdx, -hdy}};
            Eigen::MatrixXd n_0_Re{{hdy}, {hdx}};

            alpha_K[1][0][j][i] = sigma_e * (n_Ke * Lambda_K * n_1_Re)(0) / (n_0_Te * n_1_Re)(0);
            alpha_K[1][1][j][i] = sigma_e * (n_Ke * Lambda_K * n_0_Re)(0) / (n_1_Te * n_0_Re)(0);

            // West:
            double sigma_w = dy;
            Eigen::MatrixXd n_Kw{{-1, 0}}; 
            Eigen::MatrixXd n_1_Tw{{-hdx, -hdy}};
            Eigen::MatrixXd n_1_Rw{{hdy}, {-hdx}};
            Eigen::MatrixXd n_0_Tw{{-hdx, hdy}};
            Eigen::MatrixXd n_0_Rw{{-hdy}, {-hdx}};

            alpha_K[3][0][j][i] = sigma_w * (n_Kw * Lambda_K * n_1_Rw)(0) / (n_0_Tw * n_1_Rw)(0);
            alpha_K[3][1][j][i] = sigma_w * (n_Kw * Lambda_K * n_0_Rw)(0) / (n_1_Tw * n_0_Rw)(0);

            // South:
            double sigma_s = dx;
            Eigen::MatrixXd n_Ks{{0, -1}}; 
            Eigen::MatrixXd n_1_Ts{{hdx, -hdy}};
            Eigen::MatrixXd n_1_Rs{{hdy}, {hdx}};
            Eigen::MatrixXd n_0_Ts{{-hdx, -hdy}};
            Eigen::MatrixXd n_0_Rs{{hdy}, {-hdx}};

            alpha_K[2][0][j][i] = sigma_s * (n_Ks * Lambda_K * n_1_Rs)(0) / (n_0_Ts * n_1_Rs)(0);
            alpha_K[2][1][j][i] = sigma_s * (n_Ks * Lambda_K * n_0_Rs)(0) / (n_1_Ts * n_0_Rs)(0);

            // North:
            double sigma_n = dx;
            Eigen::MatrixXd n_Kn{{0, 1}}; 
            Eigen::MatrixXd n_1_Tn{{-hdx, hdy}};
            Eigen::MatrixXd n_1_Rn{{-hdy}, {-hdx}};
            Eigen::MatrixXd n_0_Tn{{hdx, hdy}};
            Eigen::MatrixXd n_0_Rn{{-hdy}, {hdx}};

            alpha_K[0][0][j][i] = sigma_n * (n_Kn * Lambda_K * n_1_Rn)(0) / (n_0_Tn * n_1_Rn)(0);
            alpha_K[0][1][j][i] = sigma_n * (n_Kn * Lambda_K * n_0_Rn)(0) / (n_1_Tn * n_0_Rn)(0);


        }
    }

    // Time loop for solving
    for (int k = 0; k < 8640; ++k) {
        // Solve using FVM solver
        // solver.solve(f, dt);
        M_ = Eigen::MatrixXd::Zero(nx*ny, nx*ny);
        S_ = Eigen::VectorXd::Zero(nx*ny);
        for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
                double p = P_MIN + hdy + dy * j;
                double a = ALPHA_LC + hdx + dx * i;
                double u1, u2, u3, u4;
                // value(u) 1, 2, 3, 4 represent North-East, North-West, South-West, South-East corner values
                if (j==0 and i==0){
                    u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u2 = 0;
                    u1 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i+1) + f(50 * (j+1) + i+1)) / 4.0;
                }
                else if (j==0 and i==nx-1){
                    u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u2 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i - 1) + f(50 * (j+1) + i-1)) / 4.0;
                    u1 = (f(50 * j + i) + f(50 * (j+1) + i)) / 2.0;
                }
                else if (j==0){
                    u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u2 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i - 1) + f(50 * (j + 1) + i - 1)) / 4.0;
                    u1 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i + 1) + f(50 * (j + 1) + i + 1)) / 4.0;
                }
                else if (i==nx-1 and j ==ny-1){
                    u1 = 0.0;
                    u2 = 0.0;
                    u3 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i - 1) + f(50 * (j - 1) + i - 1)) / 4.0;
                    u4 = (f(50 * j + i) + f(50 * (j-1) + i)) / 2.0;
                }
                else if (i==nx-1){
                    u1 = (f(50 * j + i) + f(50 * (j+1) + i)) / 2.0;
                    u4 = (f(50 * j + i) + f(50 * (j-1) + i)) / 2.0;
                    u2 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i - 1) + f(50 * (j + 1) + i - 1)) / 4.0;
                    u3 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i - 1) + f(50 * (j - 1) + i - 1)) / 4.0;
                }
                else if (j==ny-1 and i==0){
                    u1 = 0;
                    u2 = 0;
                    u3 = 0;
                    u4 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i + 1) + f(50 * (j - 1) + i + 1)) / 4.0;
                }
                else if (j=ny-1){
                    u1 = 0;
                    u2 = 0;
                    u3 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i - 1) + f(50 * (j - 1) + i - 1)) / 4.0;
                    u4 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i + 1) + f(50 * (j - 1) + i + 1)) / 4.0;
                }
                else if (i==0){
                    u2 = 0;
                    u3 = 0;
                    u1 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i + 1) + f(50 * (j + 1) + i + 1)) / 4.0;
                    u4 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i + 1) + f(50 * (j - 1) + i + 1)) / 4.0;
                }
                else{
                    u1 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i + 1) + f(50 * (j + 1) + i + 1)) / 4.0;
                    u2 = (f(50 * j + i) + f(50 * (j+1) + i) + f(50 * j + i - 1) + f(50 * (j + 1) + i - 1)) / 4.0;
                    u3 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i - 1) + f(50 * (j - 1) + i - 1)) / 4.0;
                    u4 = (f(50 * j + i) + f(50 * (j-1) + i) + f(50 * j + i + 1) + f(50 * (j - 1) + i + 1)) / 4.0;
                }

                double a_K_n = alpha_K[0][0][j][i] * u1 + alpha_K[0][1][j][i] * u2;
                double a_K_e = alpha_K[1][0][j][i] * u4 + alpha_K[1][1][j][i] * u1;
                double a_K_s = alpha_K[2][0][j][i] * u3 + alpha_K[2][1][j][i] * u4;
                double a_K_w = alpha_K[3][0][j][i] * u2 + alpha_K[3][1][j][i] * u3;

                if(j==ny-1){
                    double B_sigma = - a_K_n;
                    S_(50 * j + i) += a_K_n / G(a, p);
                    M_(50 * j + i, 50 * j + i) += (alpha_K[0][0][j][i] + alpha_K[0][1][j][i]) / G(a, p);
                }
                else{
                    double a_L_n = alpha_K[2][0][j + 1][i] * u2 + alpha_K[2][1][j + 1][i] * u1;
                    double mu_K = calMuK(a_K_n, a_L_n);
                    double mu_L = 1.0 - mu_K;
                    double B_sigma = mu_L * a_L_n - mu_K * a_K_n;
                    double B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    double B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    double A_K = mu_K * (alpha_K[0][0][j][i] + alpha_K[0][1][j][i]) + B_sigma_p / (f(50 * j + i) + 1e-15);
                    double A_L = mu_L * (alpha_K[2][0][j+1][i] + alpha_K[2][1][j+1][i]) + B_sigma_n / (f(50 * (j+1) + i) + 1e-15);
                    M_(50 * j + i, 50 * j + i) += A_K / G(a, p);
                    M_(50 * j + i, 50 * (j + 1) + i) -= A_L / G(a, p);
                }

                if(j==0){
                    double B_sigma = - a_K_s;
                    S_(50 * j + i) += a_K_s / G(a, p);
                    M_(50 * j + i, 50 * j + i) += (alpha_K[2][0][j][i] + alpha_K[2][1][j][i]) / G(a, p);
                }
                else{
                    double a_L_s = alpha_K[0][0][j-1][i] * u4 + alpha_K[0][1][j-1][i] * u3;
                    double mu_K = calMuK(a_K_s, a_L_s);
                    double mu_L = 1.0 - mu_K;
                    double B_sigma = mu_L * a_L_s - mu_K * a_K_s;
                    double B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    double B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    double A_K = mu_K * (alpha_K[2][0][j][i] + alpha_K[2][1][j][i]) + B_sigma_p / (f(50 * j + i) + 1e-15);
                    double A_L = mu_L * (alpha_K[0][0][j-1][i] + alpha_K[0][1][j-1][i]) + B_sigma_n / (f(50 * (j-1) + i) + 1e-15);
                    M_(50 * j + i, 50 * j + i) += A_K / G(a, p);
                    M_(50 * j + i, 50 * (j - 1) + i) -= A_L / G(a, p);
                }

                if (i==0){
                    double B_sigma = - a_K_w;
                    S_(50 * j + i) = a_K_w / G(a, p);
                    M_(50 * j + i, 50 * j + i) += (alpha_K[3][0][j][i] + alpha_K[3][1][j][i]) / G(a, p);
                }
                else{
                    double a_L_w = alpha_K[1][0][j][i-1] * u3 + alpha_K[1][1][j][i-1] * u2;
                    double mu_K = calMuK(a_K_w, a_L_w);
                    double mu_L = 1.0 - mu_K;
                    double B_sigma = mu_L * a_L_w - mu_K * a_K_w;
                    double B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    double B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    double A_K = mu_K * (alpha_K[3][0][j][i] + alpha_K[3][1][j][i]) + B_sigma_p / (f(50 * j + i) + 1e-15);
                    double A_L = mu_L * (alpha_K[1][0][j][i-1] + alpha_K[1][1][j][i-1]) + B_sigma_n / (f(50 * j + i - 1) + 1e-15);
                    M_(50 * j + i, 50 * j + i) += A_K / G(a, p);
                    M_(50 * j + i, 50 * j + i - 1) -= A_L / G(a, p);

                }

                if (i==nx-1){
                    continue;
                }
                else{
                    double a_L_e = alpha_K[3][0][j][i+1] * u1 + alpha_K[3][1][j][i+1] * u4;
                    double mu_K = calMuK(a_K_e, a_L_e);
                    double mu_L = 1.0 - mu_K;
                    double B_sigma = mu_L * a_L_e - mu_K * a_K_e;
                    double B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    double B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    double A_K = mu_K * (alpha_K[1][0][j][i] + alpha_K[1][1][j][i]) + B_sigma_p / (f(50 * j + i) + 1e-15);
                    double A_L = mu_L * (alpha_K[3][0][j][i+1] + alpha_K[3][1][j][i+1]) + B_sigma_n / (f(50 * j + i+1) + 1e-15);
                    M_(50 * j + i, 50 * j + i) += A_K / G(a, p);
                    M_(50 * j + i, 50 * j + i + 1) -= A_L / G(a, p);
                }
            }
        }
        Eigen::MatrixXd M_Integration = M_ * dt / (dx * dy) + Eigen::MatrixXd::Identity(nx * ny, nx * ny);
        Eigen::MatrixXd M_inv = M_Integration.inverse();
        f = M_inv * (S_ * dt / (dx * dy) + f);

        // Output or visualize f at each time step
        std::cout << "Time step: " << k << std::endl;
        
        if(k % 432 == 0){
            string path = "../output/SMPPFV_10s/smppfv" + std::to_string(int(k / 432));

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

