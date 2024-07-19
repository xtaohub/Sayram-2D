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

void FVMSolver::coeff_M_add_dirbc(int i, int j, int inbr) {  
  Edge edge; 
  m.get_nbr_edg(i, j, inbr, &edge); 

  Ind indA, indB; 
  m.indO(edge.A, &indA); 
  double fA = vertex_f_(indA.i, indA.j); 
  
  m.indO(edge.B, &indB); 
  double fB = vertex_f_(indB.i,indB.j); 

  double aK = alpha_K_(i,j,inbr).A * fA + alpha_K_(i,j,inbr).B * fB; 

  S_(ind2to1(i,j)) += aK; 
  V_coeffs_.push_back(T(ind2to1(i,j), ind2to1(i,j), (alpha_K_(i,j,inbr).A + alpha_K_(i,j,inbr).B)));
}


void FVMSolver::assemble(){ // obtain S * U^-1 and V * U^-1 
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
  
  coeff_M_add_inner(i,j,1); 
  coeff_M_add_inner(i,j,2); 
  coeff_M_add_dirbc(i, j, 0); 
  coeff_M_add_dirbc(i, j, 3); 

  // corner gird: i==nx-1, j==0
  i = nx-1;
  j = 0;

  coeff_M_add_inner(i, j, 1);
  coeff_M_add_inner(i, j, 0); 
  coeff_M_add_dirbc(i, j, 3); 

  // i min boundary, j==0
  j = 0;
  for (i=1; i<nx-1; ++i){
    coeff_M_add_inner(i, j, 1); 
    coeff_M_add_inner(i, j, 0); 
    coeff_M_add_inner(i, j, 2); 
    coeff_M_add_dirbc(i, j, 3); 
  }

  // corner gird: i==nx-1, j==ny-1
  i = nx-1;
  j = ny-1;

  //
  coeff_M_add_inner(i, j, 0);
  coeff_M_add_inner(i, j, 3);
  coeff_M_add_dirbc(i, j, 1); 

  i = nx-1;
  for(j=1; j<ny-1; j++){
    coeff_M_add_inner(i, j, 1);
    coeff_M_add_inner(i, j, 3);
    coeff_M_add_inner(i, j, 0);
  }

  // corner gird: i==0, j==ny-1
  i = 0;
  j = ny - 1;

  coeff_M_add_inner(i, j, 3);
  coeff_M_add_inner(i, j, 2);
  coeff_M_add_dirbc(i, j, 0); 
  coeff_M_add_dirbc(i, j, 1); 

  // j max boundary (except corner), j==ny-1
  j = ny - 1;
  for (i=1; i<nx-1; ++i){
    coeff_M_add_inner(i, j, 3);
    coeff_M_add_inner(i, j, 0);
    coeff_M_add_inner(i, j, 2);
    coeff_M_add_dirbc(i, j, 1); 
  }

  // i min boundary (except corner), i==0
  i = 0;
  for(j=1; j<ny-1; j++){

    coeff_M_add_inner(i, j, 1);
    coeff_M_add_inner(i, j, 3);
    coeff_M_add_inner(i, j, 2);
    coeff_M_add_dirbc(i, j, 0); 
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
