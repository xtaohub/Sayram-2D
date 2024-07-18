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
            f_(i,j) = init_f(a, p);
      }
    }
    construct_alpha_K();
}

void FVMSolver::alpha_K(const Eigen::Matrix2d& Lambda_K, const Eigen::Vector2d& K, const Eigen::Vector2d& A, const Eigen::Vector2d& B, double* alpha_KAp, double* alpha_KBp){

  Eigen::Vector2d xbk = B-K;  
  Eigen::Vector2d xak = A-K; 

  Eigen::Vector2d ba = B-A; 
  double sigma = ba.norm(); 

  Eigen::Vector2d rxbk = {xbk(1), -xbk(0)};  // (x,y) rotated by 90 clockwise, becoming (y, -x)
  Eigen::Vector2d rxak = {xak(1), -xak(0)}; 

  Eigen::Vector2d nksigma = {ba(1) / sigma, -ba(0) / sigma}; 

  *alpha_KAp = sigma * (nksigma.transpose() * Lambda_K * rxbk)(0) / (xak.transpose() * rxbk)(0);
  *alpha_KBp = sigma * (nksigma.transpose() * Lambda_K * rxak)(0) / (xbk.transpose() * rxak)(0); 
}


void FVMSolver::construct_alpha_K(){

  Eigen::Matrix2d Lambda_K;

  double a;
  double p;

  Eigen::Vector2d K;
  Eigen::Vector2d v_imjp;
  Eigen::Vector2d v_ipjp;
  Eigen::Vector2d v_imjm;
  Eigen::Vector2d v_ipjm;

  for (int i = 0; i < nx; i++){
    a = ALPHA_LC + hdx + dx * i;
    for (int j = 0; j < ny; j++){
      p = P_MIN + hdy + dy * j;

      Lambda_K << d.getDaa(0.0, i, j) * G(a, p), d.getDap(0.0, i, j) * G(a, p), 
               d.getDap(0.0, i, j) * G(a, p), d.getDpp(0.0, i, j) * G(a, p);

      K  << a, p;
      v_ipjp << a + hdx, p + hdy;
      v_imjp << a - hdx, p + hdy;
      v_imjm << a - hdx, p - hdy;
      v_ipjm << a + hdx, p - hdy;

      alpha_K(Lambda_K, K, v_ipjm, v_ipjp, &alpha_K_(i,j).ip1.A, &alpha_K_(i,j).ip1.B);
      alpha_K(Lambda_K, K, v_imjp, v_imjm, &alpha_K_(i,j).im1.A, &alpha_K_(i,j).im1.B);
      alpha_K(Lambda_K, K, v_ipjp, v_imjp, &alpha_K_(i,j).jp1.A, &alpha_K_(i,j).jp1.B);
      alpha_K(Lambda_K, K, v_imjm, v_ipjm, &alpha_K_(i,j).jm1.A, &alpha_K_(i,j).jm1.B);
    }
  }
}


void FVMSolver::coeff_M_add_n(int i, int j, double a, double p, double u_ipjp, double u_imjp){
  double a_K_n = alpha_K_(i,j).jp1.A * u_ipjp + alpha_K_(i,j).jp1.B * u_imjp;
  double a_L_n = alpha_K_(i,j+1).jm1.A * u_imjp + alpha_K_(i,j+1).jm1.B * u_ipjp;
  double mu_K = coeff_mu(a_K_n, a_L_n);
  double mu_L = 1.0 - mu_K;
  double B_sigma = mu_L * a_L_n - mu_K * a_K_n;
  double B_sigma_p = bsigma_plus(B_sigma);
  double B_sigma_n = bsigma_minus(B_sigma);
  double A_K = mu_K * (alpha_K_(i,j).jp1.A + alpha_K_(i,j).jp1.B) + B_sigma_p / (f_(i, j) + 1e-15);
  double A_L = mu_L * (alpha_K_(i,j+1).jm1.A + alpha_K_(i,j+1).jm1.B) + B_sigma_n / (f_(i, j+1) + 1e-15);
  M_coeffs_.push_back(T(nx * j + i, nx * j + i, A_K / G(a, p)));
  M_coeffs_.push_back(T(nx * j + i, nx * (j + 1) + i, -A_L / G(a, p)));
}

void FVMSolver::coeff_M_add_s(int i, int j, double a, double p, double u_imjm, double u_ipjm){
  double a_K_s = alpha_K_(i,j).jm1.A * u_imjm + alpha_K_(i,j).jm1.B * u_ipjm;
  double a_L_s = alpha_K_(i,j-1).jp1.A * u_ipjm + alpha_K_(i,j-1).jp1.B * u_imjm;
  double mu_K = coeff_mu(a_K_s, a_L_s);
  double mu_L = 1.0 - mu_K;
  double B_sigma = mu_L * a_L_s - mu_K * a_K_s;
  double B_sigma_p = bsigma_plus(B_sigma);
  double B_sigma_n = bsigma_minus(B_sigma);
  double A_K = mu_K * (alpha_K_(i,j).jm1.A + alpha_K_(i,j).jm1.B) + B_sigma_p / (f_(i, j) + 1e-15);
  double A_L = mu_L * (alpha_K_(i,j-1).jp1.A + alpha_K_(i,j-1).jp1.B) + B_sigma_n / (f_(i, j-1) + 1e-15);
  M_coeffs_.push_back(T(nx * j + i, nx * j + i, A_K / G(a, p)));
  M_coeffs_.push_back(T(nx * j + i, nx * (j - 1) + i, -A_L / G(a, p)));
}

void FVMSolver::coeff_M_add_e(int i, int j, double a, double p, double u_ipjm, double u_ipjp){
  double a_K_e = alpha_K_(i,j).ip1.A * u_ipjm + alpha_K_(i,j).ip1.B * u_ipjp;
  double a_L_e = alpha_K_(i+1, j).ip1.A * u_ipjp + alpha_K_(i+1, j).ip1.B * u_ipjm;
  double mu_K = coeff_mu(a_K_e, a_L_e);
  double mu_L = 1.0 - mu_K;
  double B_sigma = mu_L * a_L_e - mu_K * a_K_e;
  double B_sigma_p = bsigma_plus(B_sigma);
  double B_sigma_n = bsigma_minus(B_sigma);
  double A_K = mu_K * (alpha_K_(i,j).ip1.A + alpha_K_(i,j).ip1.B) + B_sigma_p / (f_(i, j) + 1e-15);
  double A_L = mu_L * (alpha_K_(i+1, j).im1.A + alpha_K_(i+1, j).im1.B) + B_sigma_n / (f_(i+1, j) + 1e-15);
  M_coeffs_.push_back(T(nx * j + i, nx * j + i, A_K / G(a, p)));
  M_coeffs_.push_back(T(nx * j + i, nx * j + i + 1, -A_L / G(a, p)));
}

void FVMSolver::coeff_M_add_w(int i, int j, double a, double p, double u_imjp, double u_imjm){
  double a_K_w = alpha_K_(i,j).im1.A * u_imjp + alpha_K_(i,j).im1.B * u_imjm;
  double a_L_w = alpha_K_(i-1, j).ip1.A * u_imjm + alpha_K_(i-1,j).ip1.B * u_imjp;
  double mu_K = coeff_mu(a_K_w, a_L_w);
  double mu_L = 1.0 - mu_K;
  double B_sigma = mu_L * a_L_w - mu_K * a_K_w;
  double B_sigma_p = bsigma_plus(B_sigma);
  double B_sigma_n = bsigma_minus(B_sigma);
  double A_K = mu_K * (alpha_K_(i,j).im1.A + alpha_K_(i,j).im1.B) + B_sigma_p / (f_(i, j) + 1e-15);
  double A_L = mu_L * (alpha_K_(i-1, j).ip1.A + alpha_K_(i-1, j).ip1.B) + B_sigma_n / (f_(i-1, j) + 1e-15);
  M_coeffs_.push_back(T(nx * j + i, nx * j + i, A_K / G(a, p)));
  M_coeffs_.push_back(T(nx * j + i, nx * j + i - 1, -A_L / G(a, p)));
}

void FVMSolver::coeff_add_n_edge(int i, int j, double a, double p, double u_ipjp, double u_imjp){
  double a_K_n = alpha_K_(i,j).jp1.A * u_ipjp + alpha_K_(i,j).jp1.B * u_imjp;
  S_(nx * j + i) += a_K_n / G(a, p);
  M_coeffs_.push_back(T(nx * j + i, nx * j + i, (alpha_K_(i,j).jp1.A + alpha_K_(i,j).jp1.B) / G(a, p)));
}

void FVMSolver::coeff_add_s_edge(int i, int j, double a, double p, double u_imjm, double u_ipjm){
  double a_K_s = alpha_K_(i,j).jm1.A * u_imjm + alpha_K_(i,j).jm1.B * u_ipjm;
  S_(nx * j + i) += a_K_s / G(a, p);
  M_coeffs_.push_back(T(nx * j + i, nx * j + i, (alpha_K_(i,j).jm1.A + alpha_K_(i,j).jm1.B) / G(a, p)));
}

void FVMSolver::coeff_add_w_edge(int i, int j, double a, double p, double u_imjp, double u_imjm){
  double a_K_w = alpha_K_(i,j).im1.A * u_imjp + alpha_K_(i,j).im1.B * u_imjm;
  S_(nx * j + i) = a_K_w / G(a, p);
  M_coeffs_.push_back(T(nx * j + i, nx * j + i, (alpha_K_(i,j).im1.A + alpha_K_(i,j).im1.B) / G(a, p)));
}


void FVMSolver::assemble(){ // obtain S and M 
  double u_ipjp, u_imjp, u_imjm, u_ipjm;
  double a, p;
  int i, j;

  // inner grids
  for (i=1; i<nx-1; ++i){
    a = ALPHA_LC + hdx + dx * i;
    for (j=1; j<ny-1; ++j){
      p = P_MIN + hdy + dy * j;
      u_ipjp = vertex_f(i,j);
      u_imjp = vertex_f(i-1,j);
      u_imjm = vertex_f(i-1, j-1);
      u_ipjm = vertex_f(i,j-1);

      coeff_M_add_n(i, j, a, p, u_ipjp, u_imjp);
      coeff_M_add_s(i, j, a, p, u_imjm, u_ipjm);
      coeff_M_add_w(i, j, a, p, u_imjp, u_imjm);
      coeff_M_add_e(i, j, a, p, u_ipjm, u_ipjp);
    }
  }

  // South-West corner，i==0, j==0, 
  i = 0;
  j = 0;
  a = ALPHA_LC + hdx + dx * i;
  p = P_MIN + hdy + dy * j;
  u_ipjp = vertex_f(i, j);
  u_imjp = 0;
  u_imjm = 0;
  u_ipjm = init_f(ALPHA_LC + dx, P_MIN);

  coeff_M_add_n(i, j, a, p, u_ipjp, u_imjp);
  coeff_M_add_e(i, j, a, p, u_ipjm, u_ipjp);
  coeff_add_w_edge(i, j, a, p, u_imjp, u_imjm);
  coeff_add_s_edge(i, j, a, p, u_imjm, u_ipjm);


  // South-East corner, i==nx-1, j==0
  i = nx-1;
  j = 0;
  a = ALPHA_LC + hdx + dx * i;
  p = P_MIN + hdy + dy * j;
  u_ipjp = (f_(i, j) + f_(i, j+1)) / 2.0;
  u_imjp = vertex_f(i-1, j);
  u_imjm = init_f(a - hdx, P_MIN);
  u_ipjm = init_f(a + hdx, P_MIN);

  coeff_M_add_n(i, j, a, p, u_ipjp, u_imjp);
  coeff_M_add_w(i, j, a, p, u_imjp, u_imjm);
  coeff_add_s_edge(i, j, a, p, u_imjm, u_ipjm);

  // South edge (except corner), j==0
  j = 0;
  p = P_MIN + hdy + dy * j;
  for (i=1; i<nx-1; ++i){
    a = ALPHA_LC + hdx + dx * i;
    u_ipjp = vertex_f(i, j);
    u_imjp = vertex_f(i-1, j);
    u_imjm = init_f(a - hdx, P_MIN);
    u_ipjm = init_f(a + hdx, P_MIN);

    coeff_M_add_n(i, j, a, p, u_ipjp, u_imjp);
    coeff_M_add_w(i, j, a, p, u_imjp, u_imjm);
    coeff_M_add_e(i, j, a, p, u_ipjm, u_ipjp);
    coeff_add_s_edge(i, j, a, p, u_imjm, u_ipjm);
  }

  // North East corner, i==nx-1, j==ny-1
  i = nx-1;
  j = ny-1;
  a = ALPHA_LC + hdx + dx * i;
  p = P_MIN + hdy + dy * j;

  u_ipjp = 0;
  u_imjp = 0;
  u_imjm = vertex_f(i-1, j-1);
  u_ipjm = (f_(i, j) + f_(i, j-1)) / 2.0;

  coeff_M_add_s(i, j, a, p, u_imjm, u_ipjm);
  coeff_M_add_w(i, j, a, p, u_imjp, u_imjm);
  coeff_add_n_edge(i, j, a, p, u_ipjp, u_imjp);

  // East edge (except corner), i==nx-1
  i = nx-1;
  a = ALPHA_LC + hdx + dx * i;
  for(j=1; j<ny-1; j++){
    p = P_MIN + hdy + dy * j;

    u_ipjp = (f_(i, j) + f_(i, j+1)) / 2.0;
    u_imjp = vertex_f(i-1, j);
    u_imjm = vertex_f(i-1, j-1);
    u_ipjm = (f_(i, j) + f_(i, j-1)) / 2.0;

    coeff_M_add_n(i, j, a, p, u_ipjp, u_imjp);
    coeff_M_add_s(i, j, a, p, u_imjm, u_ipjm);
    coeff_M_add_w(i, j, a, p, u_imjp, u_imjm);
  }

  // North West corner, i==0, j==ny-1
  i = 0;
  j = ny - 1;
  a = ALPHA_LC + hdx + dx * i;
  p = P_MIN + hdy + dy * j;

  u_ipjp = 0;
  u_imjp = 0;
  u_imjm = 0;
  u_ipjm = vertex_f(i, j-1);

  coeff_M_add_s(i, j, a, p, u_imjm, u_ipjm);
  coeff_M_add_e(i, j, a, p, u_ipjm, u_ipjp);
  coeff_add_w_edge(i, j, a, p, u_imjp, u_imjm);
  coeff_add_n_edge(i, j, a, p, u_ipjp, u_imjp);


  // North edge (except corner), j==ny-1
  j = ny - 1;
  p = P_MIN + hdy + dy * j;
  for (i=1; i<nx-1; ++i){
    a = ALPHA_LC + hdx + dx * i;

    u_ipjp = 0;
    u_imjp = 0;
    u_imjm = vertex_f(i-1, j-1);
    u_ipjm = vertex_f(i, j-1);

    coeff_M_add_s(i, j, a, p, u_imjm, u_ipjm);
    coeff_M_add_w(i, j, a, p, u_imjp, u_imjm);
    coeff_M_add_e(i, j, a, p, u_ipjm, u_ipjp);
    coeff_add_n_edge(i, j, a, p, u_ipjp, u_imjp);
  }

  // West edge (except corner), i==0
  i = 0;
  a = ALPHA_LC + hdx + dx * i;
  for(j=1; j<ny-1; j++){
    p = P_MIN + hdy + dy * j;

    u_ipjp = vertex_f(i, j);
    u_imjp = 0;
    u_imjm = 0;
    u_ipjm = vertex_f(i, j-1);

    coeff_M_add_n(i, j, a, p, u_ipjp, u_imjp);
    coeff_M_add_s(i, j, a, p, u_imjm, u_ipjm);
    coeff_M_add_e(i, j, a, p, u_ipjm, u_ipjp);
    coeff_add_w_edge(i, j, a, p, u_imjp, u_imjm);
  }

  M_.setFromTriplets(M_coeffs_.begin(), M_coeffs_.end());
}

void FVMSolver::update() {

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


