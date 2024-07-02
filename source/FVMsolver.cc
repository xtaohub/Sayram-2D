/*
 * File:        FVMsolver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */


#include "FVMSolver.hpp"
#include "Parameters.hpp"
#include "Eigen/Sparse"

typedef Eigen::Triplet<double>T;

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

void FVMSolver::initial(Eigen::VectorXd& f){
    double p;
    double a;
    for (int j = 0; j< ny; j++){
        for (int i = 0; i < nx; i++){
            p = P_MIN + hdy + dy * j;
            a = ALPHA_LC + hdx + dx * i;
            f(j * nx + i) = exp(-(p2e(p, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (p * p);
        }
    }
}

void FVMSolver::constructAlpha_K(double (&alpha_K)[4][2][ny][nx]){
    D diffusion(grid);
    diffusion.constructD(0.0);

    Eigen::Matrix2d Lambda_K;

    Eigen::Vector2d n_Ke;
    Eigen::Vector2d n_1_Te;
    Eigen::Vector2d n_1_Re;
    Eigen::Vector2d n_0_Te;
    Eigen::Vector2d n_0_Re;
    n_Ke << 1, 0; 
    n_1_Te << hdx, hdy;
    n_1_Re << -hdy, hdx;
    n_0_Te << hdx, -hdy;
    n_0_Re << hdy, hdx;

    Eigen::Vector2d n_Kw; 
    Eigen::Vector2d n_1_Tw;
    Eigen::Vector2d n_1_Rw;
    Eigen::Vector2d n_0_Tw;
    Eigen::Vector2d n_0_Rw;
    n_Kw << -1, 0;
    n_1_Tw << -hdx, -hdy;
    n_1_Rw << hdy, -hdx;
    n_0_Tw << -hdx, hdy;
    n_0_Rw << -hdy, -hdx;

                
    Eigen::Vector2d n_Ks; 
    Eigen::Vector2d n_1_Ts;
    Eigen::Vector2d n_1_Rs;
    Eigen::Vector2d n_0_Ts;
    Eigen::Vector2d n_0_Rs;
    n_Ks << 0, -1;
    n_1_Ts << hdx, -hdy;
    n_1_Rs << hdy, hdx;
    n_0_Ts << -hdx, -hdy;
    n_0_Rs << hdy, -hdx;

                
    Eigen::Vector2d n_Kn; 
    Eigen::Vector2d n_1_Tn;
    Eigen::Vector2d n_1_Rn;
    Eigen::Vector2d n_0_Tn;
    Eigen::Vector2d n_0_Rn;
    n_Kn << 0, 1; 
    n_1_Tn << -hdx, hdy;
    n_1_Rn << -hdy, -hdx;
    n_0_Tn << hdx, hdy;
    n_0_Rn << -hdy, hdx;

    double sigma_e = dy;
    double sigma_w = dy;
    double sigma_s = dx;
    double sigma_n = dx;

    for (int j = 0; j< ny; j++){
        for (int i = 0; i < nx; i++){
            // 8 direction of vector:
            double p = P_MIN + hdy + dy * j;
            double a = ALPHA_LC + hdx + dx * i;

            Lambda_K << diffusion.getDaa(0.0, j, i) * G(a, p), diffusion.getDap(0.0, j, i) * G(a, p), 
                        diffusion.getDap(0.0, j, i) * G(a, p), diffusion.getDpp(0.0, j, i) * G(a, p);
            
            // 0, 1 represent counterclockwise
            // 1 represents +1
            // e, w, n, s represent East, West, North, South edges
            
            // East:
            alpha_K[1][0][j][i] = sigma_e * (n_Ke.transpose() * Lambda_K * n_1_Re)(0) / (n_0_Te.transpose() * n_1_Re)(0);
            alpha_K[1][1][j][i] = sigma_e * (n_Ke.transpose() * Lambda_K * n_0_Re)(0) / (n_1_Te.transpose() * n_0_Re)(0);


            // West:
            alpha_K[3][0][j][i] = sigma_w * (n_Kw.transpose() * Lambda_K * n_1_Rw)(0) / (n_0_Tw.transpose() * n_1_Rw)(0);
            alpha_K[3][1][j][i] = sigma_w * (n_Kw.transpose() * Lambda_K * n_0_Rw)(0) / (n_1_Tw.transpose() * n_0_Rw)(0);

            // South:
            alpha_K[2][0][j][i] = sigma_s * (n_Ks.transpose() * Lambda_K * n_1_Rs)(0) / (n_0_Ts.transpose() * n_1_Rs)(0);
            alpha_K[2][1][j][i] = sigma_s * (n_Ks.transpose() * Lambda_K * n_0_Rs)(0) / (n_1_Ts.transpose() * n_0_Rs)(0);

            // North:
            alpha_K[0][0][j][i] = sigma_n * (n_Kn.transpose() * Lambda_K * n_1_Rn)(0) / (n_0_Tn.transpose() * n_1_Rn)(0);
            alpha_K[0][1][j][i] = sigma_n * (n_Kn.transpose() * Lambda_K * n_0_Rn)(0) / (n_1_Tn.transpose() * n_0_Rn)(0);
        }
    }
}

void FVMSolver::timeForward(Eigen::VectorXd& f, const double (&alpha_K)[4][2][ny][nx], Eigen::VectorXd &S_, std::vector<T> &M_coefficients){
    double u1, u2, u3, u4;
    double a_K_n, a_K_e, a_K_s, a_K_w;
    double mu_K, mu_L, A_K, A_L, B_sigma, B_sigma_p, B_sigma_n;
    double a_L_w, a_L_s, a_L_n, a_L_e;
    double temp0;
    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
                double p = P_MIN + hdy + dy * j;
                double a = ALPHA_LC + hdx + dx * i;
                temp0 = 0.0;
                // value(u) 1, 2, 3, 4 represent North-East, North-West, South-West, South-East corner values
                if (j==0 && i==0){
                    u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u2 = 0;
                    u1 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i+1) + f(nx * (j+1) + i+1)) / 4.0;
                }
                else if (j==0 && i==nx-1){
                    u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u2 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i - 1) + f(nx * (j+1) + i-1)) / 4.0;
                    u1 = (f(nx * j + i) + f(nx * (j+1) + i)) / 2.0;
                }
                else if (j==0){
                    u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                    u2 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i - 1) + f(nx * (j + 1) + i - 1)) / 4.0;
                    u1 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i + 1) + f(nx * (j + 1) + i + 1)) / 4.0;
                }
                else if (i==nx-1 && j ==ny-1){
                    u1 = 0.0;
                    u2 = 0.0;
                    u3 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i - 1) + f(nx * (j - 1) + i - 1)) / 4.0;
                    u4 = (f(nx * j + i) + f(nx * (j-1) + i)) / 2.0;
                }
                else if (i==nx-1){
                    u1 = (f(nx * j + i) + f(nx * (j+1) + i)) / 2.0;
                    u4 = (f(nx * j + i) + f(nx * (j-1) + i)) / 2.0;
                    u2 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i - 1) + f(nx * (j + 1) + i - 1)) / 4.0;
                    u3 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i - 1) + f(nx * (j - 1) + i - 1)) / 4.0;
                }
                else if (j==ny-1 && i==0){
                    u1 = 0.0;
                    u2 = 0.0;
                    u3 = 0.0;
                    u4 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i + 1) + f(nx * (j - 1) + i + 1)) / 4.0;
                }
                else if (j == ny-1){
                    u1 = 0.0;
                    u2 = 0.0;
                    u3 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i - 1) + f(nx * (j - 1) + i - 1)) / 4.0;
                    u4 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i + 1) + f(nx * (j - 1) + i + 1)) / 4.0;
                }
                else if (i==0){
                    u2 = 0.0;
                    u3 = 0.0;
                    u1 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i + 1) + f(nx * (j + 1) + i + 1)) / 4.0;
                    u4 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i + 1) + f(nx * (j - 1) + i + 1)) / 4.0;
                }
                else{
                    u1 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i + 1) + f(nx * (j + 1) + i + 1)) / 4.0;
                    u2 = (f(nx * j + i) + f(nx * (j+1) + i) + f(nx * j + i - 1) + f(nx * (j + 1) + i - 1)) / 4.0;
                    u3 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i - 1) + f(nx * (j - 1) + i - 1)) / 4.0;
                    u4 = (f(nx * j + i) + f(nx * (j-1) + i) + f(nx * j + i + 1) + f(nx * (j - 1) + i + 1)) / 4.0;
                }

                a_K_n = alpha_K[0][0][j][i] * u1 + alpha_K[0][1][j][i] * u2;
                a_K_e = alpha_K[1][0][j][i] * u4 + alpha_K[1][1][j][i] * u1;
                a_K_s = alpha_K[2][0][j][i] * u3 + alpha_K[2][1][j][i] * u4;
                a_K_w = alpha_K[3][0][j][i] * u2 + alpha_K[3][1][j][i] * u3;

                if(j==ny-1){
                    B_sigma = - a_K_n;
                    S_(nx * j + i) += a_K_n / G(a, p);
                    temp0 += (alpha_K[0][0][j][i] + alpha_K[0][1][j][i]) / G(a, p);
                }
                else{
                    a_L_n = alpha_K[2][0][j + 1][i] * u2 + alpha_K[2][1][j + 1][i] * u1;
                    mu_K = calMuK(a_K_n, a_L_n);
                    mu_L = 1.0 - mu_K;
                    B_sigma = mu_L * a_L_n - mu_K * a_K_n;
                    B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    A_K = mu_K * (alpha_K[0][0][j][i] + alpha_K[0][1][j][i]) + B_sigma_p / (f(nx * j + i) + 1e-15);
                    A_L = mu_L * (alpha_K[2][0][j+1][i] + alpha_K[2][1][j+1][i]) + B_sigma_n / (f(nx * (j+1) + i) + 1e-15);
                    temp0 += A_K / G(a, p);
                    M_coefficients.push_back(T(nx * j + i, nx * (j + 1) + i, -A_L / G(a, p)));
                }

                if(j==0){
                    B_sigma = - a_K_s;
                    S_(nx * j + i) += a_K_s / G(a, p);
                    temp0 += (alpha_K[2][0][j][i] + alpha_K[2][1][j][i]) / G(a, p);
                }
                else{
                    a_L_s = alpha_K[0][0][j-1][i] * u4 + alpha_K[0][1][j-1][i] * u3;
                    mu_K = calMuK(a_K_s, a_L_s);
                    mu_L = 1.0 - mu_K;
                    B_sigma = mu_L * a_L_s - mu_K * a_K_s;
                    B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    A_K = mu_K * (alpha_K[2][0][j][i] + alpha_K[2][1][j][i]) + B_sigma_p / (f(nx * j + i) + 1e-15);
                    A_L = mu_L * (alpha_K[0][0][j-1][i] + alpha_K[0][1][j-1][i]) + B_sigma_n / (f(nx * (j-1) + i) + 1e-15);
                    temp0 += A_K / G(a, p);
                    M_coefficients.push_back(T(nx * j + i, nx * (j - 1) + i, -A_L / G(a, p)));
                }

                if (i==0){
                    B_sigma = - a_K_w;
                    S_(nx * j + i) = a_K_w / G(a, p);
                    temp0 += (alpha_K[3][0][j][i] + alpha_K[3][1][j][i]) / G(a, p);   
                }
                else{
                    a_L_w = alpha_K[1][0][j][i-1] * u3 + alpha_K[1][1][j][i-1] * u2;
                    mu_K = calMuK(a_K_w, a_L_w);
                    mu_L = 1.0 - mu_K;
                    B_sigma = mu_L * a_L_w - mu_K * a_K_w;
                    B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    A_K = mu_K * (alpha_K[3][0][j][i] + alpha_K[3][1][j][i]) + B_sigma_p / (f(nx * j + i) + 1e-15);
                    A_L = mu_L * (alpha_K[1][0][j][i-1] + alpha_K[1][1][j][i-1]) + B_sigma_n / (f(nx * j + i - 1) + 1e-15);
                    temp0 += A_K / G(a, p);
                    M_coefficients.push_back(T(nx * j + i, nx * j + i - 1, -A_L / G(a, p)));
                }

                if (i==nx-1){
                    M_coefficients.push_back(T(nx * j + i, nx * j + i, temp0));
                    continue;
                }
                else{
                    a_L_e = alpha_K[3][0][j][i+1] * u1 + alpha_K[3][1][j][i+1] * u4;
                    mu_K = calMuK(a_K_e, a_L_e);
                    mu_L = 1.0 - mu_K;
                    B_sigma = mu_L * a_L_e - mu_K * a_K_e;
                    B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                    B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                    A_K = mu_K * (alpha_K[1][0][j][i] + alpha_K[1][1][j][i]) + B_sigma_p / (f(nx * j + i) + 1e-15);
                    A_L = mu_L * (alpha_K[3][0][j][i+1] + alpha_K[3][1][j][i+1]) + B_sigma_n / (f(nx * j + i+1) + 1e-15);
                    temp0 += A_K / G(a, p);
                    M_coefficients.push_back(T(nx * j + i, nx * j + i, temp0));
                    M_coefficients.push_back(T(nx * j + i, nx * j + i + 1, -A_L / G(a, p)));
                }
            }
        }
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


