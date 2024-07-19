/*
 * File:        FVMsolver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */


#include "FVMSolver.hpp"

FVMSolver::FVMSolver(const Mesh& m_in, const D& d_in, const BoundaryConditions& bc_in)
  : m(m_in), d(d_in), bc(bc_in){

    int nx = m.nx();
    int ny = m.ny();

    V_.resize(nx*ny, nx*ny);
    M_.resize(nx*ny, nx*ny);
    f_.resize(nx, ny);
    S_.resize(nx*ny); 
    R_.resize(nx*ny);

    alpha_K_.resize(nx, ny, m.nnbrs()); 
    vertex_f_.resize(nx+1,ny+1); 
    U_.resize(nx*ny, nx*ny);

    hdx_ = m.dx()/2.0;
    hdy_ = m.dy()/2.0;

    initial(); 

  }

void FVMSolver::initial(){
  double a;
  double p;
  for (int i = 0; i < nx; i++){
    a = m.x(i);
    for (int j = 0; j < ny; j++){
      p = m.y(j);
      f_(i,j) = init_f(a, p);
    }
  }
  construct_alpha_K();

  set_vertex_f();

  construct_U_();
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

  double a, p;

  Eigen::Vector2d K;

  Edge edge;  

  for (int i = 0; i < nx; i++){
    a = m.x(i);
    for (int j = 0; j < ny; j++){
      p = m.y(j);

      Lambda_K << d.getDaa(0.0, i, j) * G(a, p), d.getDap(0.0, i, j) * G(a, p), 
               d.getDap(0.0, i, j) * G(a, p), d.getDpp(0.0, i, j) * G(a, p);

      K  << a, p;

      for (int inbr=0; inbr<m.nnbrs(); ++inbr){
        m.get_nbr_edg(i, j, inbr, &edge);  
        alpha_K(Lambda_K, K, edge.A, edge.B, &alpha_K_(i,j, inbr).A, &alpha_K_(i,j, inbr).B); 
      }
    }
  }
}

void FVMSolver::coeff_M_add_inner(int i, int j, int inbr){
  Ind ind; 
  Edge edge; 
  m.get_nbr_edg(i, j, inbr, &edge); 
  m.get_nbr_ind(i, j, inbr, &ind); 

  Ind indA, indB; 
  m.indO(edge.A, &indA); 
  double fA = vertex_f_(indA.i, indA.j); 
  
  m.indO(edge.B, &indB); 
  double fB = vertex_f_(indB.i,indB.j); 

  double aK = coeff_a(alpha_K_(i,j,inbr).A, fA, alpha_K_(i,j,inbr).B, fB);

  int rinbr = m.rinbr(inbr); // for the neighboring cell, A and B are reversed. 
  double aL = coeff_a(alpha_K_(ind.i,ind.j,rinbr).A, fB, alpha_K_(ind.i,ind.j,rinbr).B, fA); 

  double muK = coeff_mu(aK, aL); 
  double muL = 1.0 - muK; 

  double B_sigma = muL * aL - muK * aK;
  double B_sigma_p = bsigma_plus(B_sigma);
  double B_sigma_n = bsigma_minus(B_sigma);

  double A_K = muK * (alpha_K_(i,j, inbr).A + alpha_K_(i,j, inbr).B) + B_sigma_p / (f_(i, j) + 1e-15);
  double A_L = muL * (alpha_K_(ind.i,ind.j, rinbr).A + alpha_K_(ind.i,ind.j, inbr).B) + B_sigma_n / (f_(ind.i, ind.j) + 1e-15);

  V_coeffs_.push_back(T(ind2to1(i,j), ind2to1(i,j), A_K));
  V_coeffs_.push_back(T(ind2to1(i,j), ind2to1(ind.i,ind.j), -A_L));

}

void FVMSolver::coeff_add_jp_edge(int i, int j, double u_ipjp, double u_imjp){
  double a_K = alpha_K_(i,j, 1).A * u_ipjp + alpha_K_(i,j,1).B * u_imjp;
  S_(ind2to1(i,j)) += a_K;
  V_coeffs_.push_back(T(ind2to1(i,j), ind2to1(i,j), (alpha_K_(i,j, 1).A + alpha_K_(i,j,1).B)));
}

void FVMSolver::coeff_add_jm_edge(int i, int j, double u_imjm, double u_ipjm){
  double a_K = alpha_K_(i,j,3).A * u_imjm + alpha_K_(i,j,3).B * u_ipjm;
  S_(ind2to1(i,j)) += a_K;
  V_coeffs_.push_back(T(ind2to1(i,j), ind2to1(i,j), (alpha_K_(i,j,3).A + alpha_K_(i,j,3).B)));
}

void FVMSolver::coeff_add_im_edge(int i, int j, double u_imjp, double u_imjm){
  double a_K = alpha_K_(i,j,0).A * u_imjp + alpha_K_(i,j,0).B * u_imjm;
  S_(ind2to1(i,j)) = a_K;
  V_coeffs_.push_back(T(ind2to1(i,j), ind2to1(i,j), (alpha_K_(i,j,0).A + alpha_K_(i,j,0).B)));
}


void FVMSolver::assemble(){ // obtain S * U^-1 and V * U^-1 
  double u_ipjp, u_imjp, u_imjm, u_ipjm;
  double a, p;
  int i, j;

  // inner grids
  for (i=1; i<m.nx()-1; ++i){
    for (j=1; j<m.ny()-1; ++j){

      for (int inbr = 0; inbr < m.nnbrs(); ++inbr) 
        coeff_M_add_inner(i, j, inbr);
    }
  }

  // corner gird: i==0, j==0, 
  i = 0;
  j = 0;
  a = m.x(i);
  p = m.y(j);
  u_ipjp = vertex_f_(1, 1);
  u_imjp = boundary_imin(p + hdy);
  u_imjm = boundary_imin(p + hdy);
  u_ipjm = boundary_jmin(a + hdx);
  
  coeff_M_add_inner(i,j,1); 
  coeff_M_add_inner(i,j,2); 
  coeff_add_im_edge(i, j, u_imjp, u_imjm);
  coeff_add_jm_edge(i, j, u_imjm, u_ipjm);


  // corner gird: i==nx-1, j==0
  i = nx-1;
  j = 0;
  a = m.x(i);
  p = m.y(j);
  u_ipjp = (f_(i, j) + f_(i, j+1)) / 2.0;
  u_imjp = vertex_f_(i, j+1);
  u_imjm = boundary_jmin(a - hdx);
  u_ipjm = boundary_jmin(a + hdx);

  coeff_M_add_inner(i, j, 1);
  coeff_M_add_inner(i, j, 0); 
  coeff_add_jm_edge(i, j, u_imjm, u_ipjm);

  // i min boundary, j==0
  j = 0;
  p = m.y(j);
  for (i=1; i<nx-1; ++i){
    a = m.x(i);
    u_ipjp = vertex_f_(i+1, j+1);
    u_imjp = vertex_f_(i, j+1);
    u_imjm = boundary_jmin(a - hdx);
    u_ipjm = boundary_jmin(a + hdx);

    coeff_M_add_inner(i, j, 1); 
    coeff_M_add_inner(i, j, 0); 
    coeff_M_add_inner(i, j, 2); 
    coeff_add_jm_edge(i, j, u_imjm, u_ipjm);
  }

  // corner gird: i==nx-1, j==ny-1
  i = nx-1;
  j = ny-1;
  a = m.x(i);
  p = m.y(j);

  u_ipjp = boundary_jmax(a + hdx);
  u_imjp = boundary_jmax(a - hdx);
  u_imjm = vertex_f_(i, j);
  u_ipjm = (f_(i, j) + f_(i, j-1)) / 2.0;

  //
  coeff_M_add_inner(i, j, 0);
  coeff_M_add_inner(i, j, 3);
  coeff_add_jp_edge(i, j, u_ipjp, u_imjp);

  // i max boundary (except corner), i==nx-1
  i = nx-1;
  a = m.x(i);
  for(j=1; j<ny-1; j++){
    p = m.y(j);

    u_ipjp = (f_(i, j) + f_(i, j+1)) / 2.0;
    u_imjp = vertex_f_(i, j+1);
    u_imjm = vertex_f_(i, j);
    u_ipjm = (f_(i, j) + f_(i, j-1)) / 2.0;

    coeff_M_add_inner(i, j, 1);
    coeff_M_add_inner(i, j, 3);
    coeff_M_add_inner(i, j, 0);
  }

  // corner gird: i==0, j==ny-1
  i = 0;
  j = ny - 1;
  a = m.x(i);
  p = m.y(j);

  u_ipjp = boundary_jmax(a + hdx);
  u_imjp = boundary_jmax(a - hdx);
  u_imjm = boundary_imin(p - hdy);
  u_ipjm = vertex_f_(i+1, j);

  // coeff_add_jm(i, j, a, p, u_imjm, u_ipjm);
  // coeff_add_ip(i, j, a, p, u_ipjm, u_ipjp);
  coeff_M_add_inner(i, j, 3);
  coeff_M_add_inner(i, j, 2);

  coeff_add_im_edge(i, j, u_imjp, u_imjm);
  coeff_add_jp_edge(i, j, u_ipjp, u_imjp);


  // j max boundary (except corner), j==ny-1
  j = ny - 1;
  p = m.y(j);
  for (i=1; i<nx-1; ++i){
    a = m.x(i);

    u_ipjp = boundary_jmax(a + hdx);
    u_imjp = boundary_jmax(a - hdx);
    u_imjm = vertex_f(i-1, j-1);
    u_ipjm = vertex_f(i, j-1);

    coeff_M_add_inner(i, j, 3);
    coeff_M_add_inner(i, j, 0);
    coeff_M_add_inner(i, j, 2);
    coeff_add_jp_edge(i, j, u_ipjp, u_imjp);
  }

  // i min boundary (except corner), i==0
  i = 0;
  a = m.x(i);
  for(j=1; j<ny-1; j++){
    p = m.y(j);

    u_ipjp = vertex_f(i, j);
    u_imjp = boundary_imin(p + hdy);
    u_imjm = boundary_imin(p - hdy);
    u_ipjm = vertex_f(i, j-1);

    coeff_M_add_inner(i, j, 1);
    coeff_M_add_inner(i, j, 3);
    coeff_M_add_inner(i, j, 2);

    coeff_add_im_edge(i, j, u_imjp, u_imjm);
  }

  V_.setFromTriplets(V_coeffs_.begin(), V_coeffs_.end());
}

void FVMSolver::construct_U_(){
  double a, p;
  for (int i=0; i<m.nx(); ++i){
    a = m.x(i);
    for (int j=0; j<m.ny(); ++j){
      p = m.y(j);
      U_(ind2to1(i, j), ind2to1(i, j)) = G(a, p) * m.area_dt();
    }
  }
}


void FVMSolver::update() {

  S_.setZero();
  V_coeffs_.clear();

  assemble(); 

  M_ = V_ + U_;
  R_ = S_ + U_ * f_.reshaped();

  solver.analyzePattern(M_);
  solver.factorize(M_);

  f_.reshaped() = solver.solve(R_);

  set_vertex_f();
}

void FVMSolver::set_vertex_f(){
  for (int i=1; i<m.nx(); ++i)
    for (int j=1; j<m.ny(); ++j){
      vertex_f_(i,j) = (f_(i-1,j-1) + f_(i-1,j) + f_(i,j-1) + f_(i,j)) / 4.0; 
    }

  // i == 0 boundary
  double a, p;
  for (int j = 0; j<vertex_f_.ncols(); ++j){
    p = m.yO()+j*m.dy(); 
    vertex_f_(0, j) = boundary_imin(p);
  }

  // j == 0  and j == m.ny() boundary
  for (int i = 0; i < vertex_f_.nrows(); ++i) {
    a = m.xO() + i*m.dx(); 
    vertex_f_(i,0) = boundary_jmin(a);
    vertex_f_(i,m.ny()) = boundary_jmax(a);  
  }

}
