/*
 * File:        FVMsolver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */


// #include "mkl.h"
#include "FVMSolver.hpp"

FVMSolver::FVMSolver(const Grid& g_in, const D& d_in, const BoundaryConditions& bc_in)
    : g(g_in), d(d_in), bc(bc_in){

    int nx = g.nx();
    int ny = g.ny();

    M_.resize(nx*ny, nx*ny);
    f_.resize(nx, ny);
    S_.resize(nx*ny); 

    alpha_K_.resize(nx, ny); 

    Id_ = Eigen::MatrixXd::Identity(nx*ny, nx*ny);

    hdx_ = g.dx()/2.0;
    hdy_ = g.dy()/2.0;

}

void FVMSolver::initial(){
    double a;
    double p;
    for (int i = 0; i < nx; i++){
        a = ALPHA_LC + hdx + dx * i;
        for (int j = 0; j < ny; j++){
            p = P_MIN + hdy + dy * j;
            f_(i,j) = exp(-(p2e(p, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (p * p);
      }
    }
    construct_alpha_K();
}

void FVMSolver::construct_alpha_K(){

    Eigen::Matrix2d Lambda_K;

    Eigen::Vector2d n_Ke;
    Eigen::Vector2d n_1_Te;
    Eigen::Vector2d n_1_Re;
    Eigen::Vector2d n_0_Te;
    Eigen::Vector2d n_0_Re;
    n_Ke << 1, 0; 
    n_1_Te << hdx, hdy;
    n_1_Re << hdy, -hdx;
    n_0_Te << hdx, -hdy;
    n_0_Re << -hdy, -hdx;

    Eigen::Vector2d n_Kw; 
    Eigen::Vector2d n_1_Tw;
    Eigen::Vector2d n_1_Rw;
    Eigen::Vector2d n_0_Tw;
    Eigen::Vector2d n_0_Rw;
    n_Kw << -1, 0;
    n_1_Tw << -hdx, -hdy;
    n_1_Rw << -hdy, hdx;
    n_0_Tw << -hdx, hdy;
    n_0_Rw << hdy, hdx;

    Eigen::Vector2d n_Kn; 
    Eigen::Vector2d n_1_Tn;
    Eigen::Vector2d n_1_Rn;
    Eigen::Vector2d n_0_Tn;
    Eigen::Vector2d n_0_Rn;
    n_Kn << 0, 1; 
    n_1_Tn << -hdx, hdy;
    n_1_Rn << hdy, hdx;
    n_0_Tn << hdx, hdy;
    n_0_Rn << hdy, -hdx;

    Eigen::Vector2d n_Ks; 
    Eigen::Vector2d n_1_Ts;
    Eigen::Vector2d n_1_Rs;
    Eigen::Vector2d n_0_Ts;
    Eigen::Vector2d n_0_Rs;
    n_Ks << 0, -1;
    n_1_Ts << hdx, -hdy;
    n_1_Rs << -hdy, -hdx;
    n_0_Ts << -hdx, -hdy;
    n_0_Rs << -hdy, hdx;

    double sigma_e = dy;
    double sigma_w = dy;
    double sigma_n = dx;
    double sigma_s = dx;

    double a;
    double p;

    for (int i = 0; i < nx; i++){
        a = ALPHA_LC + hdx + dx * i;
        for (int j = 0; j < ny; j++){
            p = P_MIN + hdy + dy * j;
            Lambda_K << d.getDaa(0.0, i, j) * G(a, p), d.getDap(0.0, i, j) * G(a, p), 
                        d.getDap(0.0, i, j) * G(a, p), d.getDpp(0.0, i, j) * G(a, p);
            
            // 0, 1 are counterclockwise vertex marker
            // T: Transpose vector
            // R: Rotate a vector 90 degree clockwise
            // the north represent the growth direction of p, and east is the alpha growing direction

            // East:
            alpha_K_(i,j).e.A = sigma_e * (n_Ke.transpose() * Lambda_K * n_1_Re)(0) / (n_0_Te.transpose() * n_1_Re)(0);
            alpha_K_(i,j).e.B = sigma_e * (n_Ke.transpose() * Lambda_K * n_0_Re)(0) / (n_1_Te.transpose() * n_0_Re)(0);

            // West:
            alpha_K_(i,j).w.A = sigma_w * (n_Kw.transpose() * Lambda_K * n_1_Rw)(0) / (n_0_Tw.transpose() * n_1_Rw)(0);
            alpha_K_(i,j).w.B = sigma_w * (n_Kw.transpose() * Lambda_K * n_0_Rw)(0) / (n_1_Tw.transpose() * n_0_Rw)(0);

            // North:
            alpha_K_(i,j).n.A = sigma_n * (n_Kn.transpose() * Lambda_K * n_1_Rn)(0) / (n_0_Tn.transpose() * n_1_Rn)(0);
            alpha_K_(i,j).n.B = sigma_n * (n_Kn.transpose() * Lambda_K * n_0_Rn)(0) / (n_1_Tn.transpose() * n_0_Rn)(0);

            // South:
            alpha_K_(i,j).s.A = sigma_s * (n_Ks.transpose() * Lambda_K * n_1_Rs)(0) / (n_0_Ts.transpose() * n_1_Rs)(0);
            alpha_K_(i,j).s.B = sigma_s * (n_Ks.transpose() * Lambda_K * n_0_Rs)(0) / (n_1_Ts.transpose() * n_0_Rs)(0);
        }
    }
}

void FVMSolver::assemble(){ // obtain S and M 
    double u1, u2, u3, u4;
    double a_K_n, a_K_e, a_K_s, a_K_w;
    double mu_K, mu_L, A_K, A_L, B_sigma, B_sigma_p, B_sigma_n;
    double a_L_w, a_L_s, a_L_n, a_L_e;
    double temp0;
    double a;
    double p;
    for (int i = 0; i < nx; i++){
        a = ALPHA_LC + hdx + dx * i;
        for (int j = 0; j < ny; j++){
            p = P_MIN + hdy + dy * j;
            temp0 = 0.0;
            // value(u) 1, 2, 3, 4 represent North-East, North-West, South-West, South-East corner values
            // cases when grid on corner, edge or inner the domain
            if (i==0 && j==0){
                u1 = (f_(i, j) + f_(i, j+1) + f_(i+1, j) + f_(i+1, j+1)) / 4.0;
                u2 = 0;
                u3 = 0;
                u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
            }
            else if (i==nx-1 && j==0){
                u1 = (f_(i, j) + f_(i, j+1)) / 2.0;
                u2 = (f_(i, j) + f_(i, j+1) + f_(i-1, j) + f_(i-1, j+1)) / 4.0;
                u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
            }
            else if (j==0){
                u1 = (f_(i, j) + f_(i, j+1) + f_(i+1, j) + f_(i+1, j+1)) / 4.0;
                u2 = (f_(i, j) + f_(i, j+1) + f_(i-1, j) + f_(i-1, j+1)) / 4.0;
                u3 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a - hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
                u4 = exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a + hdx) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
            }
            else if (i==nx-1 && j==ny-1){
                u1 = 0;
                u2 = 0;
                u3 = (f_(i, j) + f_(i, j-1) + f_(i-1, j) + f_(i-1, j-1)) / 4.0;
                u4 = (f_(i, j) + f_(i, j-1)) / 2.0;
            }
            else if (i==nx-1){
                u1 = (f_(i, j) + f_(i, j+1)) / 2.0;
                u2 = (f_(i, j) + f_(i, j+1) + f_(i-1, j) + f_(i-1, j+1)) / 4.0;
                u3 = (f_(i, j) + f_(i, j-1) + f_(i-1, j) + f_(i-1, j-1)) / 4.0;
                u4 = (f_(i, j) + f_(i, j-1)) / 2.0;
            }
            else if (i==0 && j==ny-1){
                u1 = 0;
                u2 = 0;
                u3 = 0;
                u4 = (f_(i, j) + f_(i, j-1) + f_(i+1, j) + f_(i+1, j-1)) / 4.0;
            }
            else if (j==ny-1){
                u1 = 0;
                u2 = 0;
                u3 = (f_(i, j) + f_(i, j-1) + f_(i-1, j) + f_(i-1, j-1)) / 4.0;
                u4 = (f_(i, j) + f_(i, j-1) + f_(i+1, j) + f_(i+1, j-1)) / 4.0;
            }
            else if (i==0){
                u1 = (f_(i, j) + f_(i, j+1) + f_(i+1, j) + f_(i+1, j+1)) / 4.0;
                u2 = 0;
                u3 = 0;
                u4 = (f_(i, j) + f_(i, j-1) + f_(i+1, j) + f_(i+1, j-1)) / 4.0;
            }
            else{
                u1 = (f_(i, j) + f_(i, j+1) + f_(i+1, j) + f_(i+1, j+1)) / 4.0;
                u2 = (f_(i, j) + f_(i, j+1) + f_(i-1, j) + f_(i-1, j+1)) / 4.0;
                u3 = (f_(i, j) + f_(i, j-1) + f_(i-1, j) + f_(i-1, j-1)) / 4.0;
                u4 = (f_(i, j) + f_(i, j-1) + f_(i+1, j) + f_(i+1, j-1)) / 4.0;
            }

            a_K_n = alpha_K_(i,j).n.A * u1 + alpha_K_(i,j).n.B * u2;
            a_K_w = alpha_K_(i,j).w.A * u2 + alpha_K_(i,j).w.B * u3;
            a_K_s = alpha_K_(i,j).s.A * u3 + alpha_K_(i,j).s.B * u4;
            a_K_e = alpha_K_(i,j).e.A * u4 + alpha_K_(i,j).e.B * u1;

            if(j==ny-1){
                B_sigma = - a_K_n;
                S_(nx * j + i) += a_K_n / G(a, p);
                temp0 += (alpha_K_(i,j).n.A + alpha_K_(i,j).n.B) / G(a, p);
            }
            else{
                a_L_n = alpha_K_(i,j+1).s.A * u2 + alpha_K_(i,j+1).s.B * u1;
                mu_K = calMuK(a_K_n, a_L_n);
                mu_L = 1.0 - mu_K;
                B_sigma = mu_L * a_L_n - mu_K * a_K_n;
                B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                A_K = mu_K * (alpha_K_(i,j).n.A + alpha_K_(i,j).n.B) + B_sigma_p / (f_(i, j) + 1e-15);
                A_L = mu_L * (alpha_K_(i,j+1).s.A + alpha_K_(i,j+1).s.B) + B_sigma_n / (f_(i, j+1) + 1e-15);
                temp0 += A_K / G(a, p);
                M_coeffs_.push_back(T(nx * j + i, nx * (j + 1) + i, -A_L / G(a, p)));
            }
            
            if(j==0){
                B_sigma = - a_K_s;
                S_(nx * j + i) += a_K_s / G(a, p);
                temp0 += (alpha_K_(i,j).s.A + alpha_K_(i,j).s.B) / G(a, p);
            }
            else{
                a_L_s = alpha_K_(i,j-1).n.A * u4 + alpha_K_(i,j-1).n.B * u3;
                mu_K = calMuK(a_K_s, a_L_s);
                mu_L = 1.0 - mu_K;
                B_sigma = mu_L * a_L_s - mu_K * a_K_s;
                B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                A_K = mu_K * (alpha_K_(i,j).s.A + alpha_K_(i,j).s.B) + B_sigma_p / (f_(i, j) + 1e-15);
                A_L = mu_L * (alpha_K_(i,j-1).n.A + alpha_K_(i,j-1).n.B) + B_sigma_n / (f_(i, j-1) + 1e-15);
                temp0 += A_K / G(a, p);
                M_coeffs_.push_back(T(nx * j + i, nx * (j - 1) + i, -A_L / G(a, p)));
            }
            
            if (i==0){
                B_sigma = - a_K_w;
                S_(nx * j + i) = a_K_w / G(a, p);
                temp0 += (alpha_K_(i,j).w.A + alpha_K_(i,j).w.B) / G(a, p);   
            }
            else{
                a_L_w = alpha_K_(i-1, j).e.A * u3 + alpha_K_(i-1,j).e.B * u2;
                mu_K = calMuK(a_K_w, a_L_w);
                mu_L = 1.0 - mu_K;
                B_sigma = mu_L * a_L_w - mu_K * a_K_w;
                B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                A_K = mu_K * (alpha_K_(i,j).w.A + alpha_K_(i,j).w.B) + B_sigma_p / (f_(i, j) + 1e-15);
                A_L = mu_L * (alpha_K_(i-1, j).e.A + alpha_K_(i-1, j).e.B) + B_sigma_n / (f_(i-1, j) + 1e-15);
                temp0 += A_K / G(a, p);
                M_coeffs_.push_back(T(nx * j + i, nx * j + i - 1, -A_L / G(a, p)));
            }
            
            if (i==nx-1){
                M_coeffs_.push_back(T(nx * j + i, nx * j + i, temp0));
                continue; // Neumann boundary condition
            }
            else{
                a_L_e = alpha_K_(i+1, j).e.A * u1 + alpha_K_(i+1, j).e.B * u4;
                mu_K = calMuK(a_K_e, a_L_e);
                mu_L = 1.0 - mu_K;
                B_sigma = mu_L * a_L_e - mu_K * a_K_e;
                B_sigma_p = (abs(B_sigma) + B_sigma) / 2;
                B_sigma_n = (abs(B_sigma) - B_sigma) / 2;
                A_K = mu_K * (alpha_K_(i,j).e.A + alpha_K_(i,j).e.B) + B_sigma_p / (f_(i, j) + 1e-15);
                A_L = mu_L * (alpha_K_(i+1, j).w.A + alpha_K_(i+1, j).w.B) + B_sigma_n / (f_(i+1, j) + 1e-15);
                temp0 += A_K / G(a, p);
                M_coeffs_.push_back(T(nx * j + i, nx * j + i, temp0));
                M_coeffs_.push_back(T(nx * j + i, nx * j + i + 1, -A_L / G(a, p)));
            }
        }
    }

    M_.setFromTriplets(M_coeffs_.begin(), M_coeffs_.end());
}

void FVMSolver::update() {
    // assembleSystem(f, d, dt);

    S_.setZero();
    M_coeffs_.clear();

    assemble(); 
 
    double dt_area = g.dt() / (g.dx() * g.dy());
    M_ = M_ * dt_area + Id_;
    S_ = S_ * dt_area + f_.reshaped();

    solver.analyzePattern(M_);
    solver.factorize(M_);

    f_.reshaped() = solver.solve(S_);
}


