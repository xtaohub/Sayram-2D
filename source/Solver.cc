/*
 * File:        Solver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */


#include "Solver.h"
#include "Parameters.h"
#include <limits>

Solver::Solver(const Parameters& paras_in, const Mesh& m_in, const D& d_in, const BCs& bcs_in)
  : paras(paras_in), m(m_in), d(d_in), bcs(bcs_in){

    std::size_t nx = m.nx();
    std::size_t ny = m.ny();

    M_.resize(nx*ny,nx*ny);

    f_.resize(nx,ny);
    R_.resize(nx*ny);

    alpha_osf_.resize({nx, ny, m.nnbrs()});
    vertex_f_.resize({nx+1, ny+1});

    init();
  }

void Solver::init(){
  double a0;
  double p;

  for (std::size_t i = 0; i < m.nx(); i++){
    a0 = m.x(i);
    for (std::size_t j = 0; j < m.ny(); j++){
      p = m.y(j);
      f_(i,j) = bcs.init_f(a0, p);
    }
  }

  t_ = 0;
  construct_alpha_osf();
  update_vertex_f();
}

void Solver::alpha_osf_func(const Eigen::Matrix2d& Lambda_K, const Point& K, const Point& A, const Point& B, NTPFA_node* nodep){

  Eigen::Vector2d xbk = B-K;  
  Eigen::Vector2d xak = A-K; 

  Eigen::Vector2d ba = B-A; 
  double sigma = ba.norm(); 

  Eigen::Vector2d rxbk = {xbk(1), -xbk(0)};  // (x,y) rotated by 90 clockwise, becoming (y, -x)
  Eigen::Vector2d rxak = {xak(1), -xak(0)}; 

  Eigen::Vector2d nksigma = {ba(1) / sigma, -ba(0) / sigma}; 

  nodep->A = sigma * (nksigma.transpose() * Lambda_K * rxbk)(0) / (xak.transpose() * rxbk)(0);
  nodep->B = sigma * (nksigma.transpose() * Lambda_K * rxak)(0) / (xbk.transpose() * rxak)(0); 
}


void Solver::construct_alpha_osf(){

  Eigen::Matrix2d Lambda_K;

  double a0, p;
  Point K;
  Edge edge;  

  for (std::size_t i = 0; i < m.nx(); i++){
    a0 = m.x(i);
    for (std::size_t j = 0; j < m.ny(); j++){
      p = m.y(j);

      Lambda_K << d.Daa(t(), i, j) * G(a0, p), d.Dap(t(), i, j) * G(a0, p), 
               d.Dap(t(), i, j) * G(a0, p), d.Dpp(t(), i, j) * G(a0, p);

      K  << a0, p;

      for (size_t inbr=0; inbr<m.nnbrs(); ++inbr){
        m.get_nbr_edg(i, j, inbr, &edge);  
        alpha_osf_func(Lambda_K, K, edge.A, edge.B, &alpha_osf_(i,j, inbr)); 
      }
    }
  }
}

void Solver::coeff_add_inner(int i, int j, int inbr){
  Ind ind; 
  Edge edge; 
  m.get_nbr_edg(i, j, inbr, &edge); 
  m.get_nbr_ind(i, j, inbr, &ind); 

  Ind indA, indB; 
  m.indO(edge.A, &indA); 
  double fA = vertex_f_(indA.i, indA.j); 

  m.indO(edge.B, &indB); 
  double fB = vertex_f_(indB.i,indB.j); 

  double aK = alpha_osf_(i,j,inbr).A * fA + alpha_osf_(i,j,inbr).B * fB;

  int rinbr = m.rinbr(inbr); // for the neighboring cell, A and B are reversed. 
  double aL = alpha_osf_(ind.i,ind.j,rinbr).A * fB + alpha_osf_(ind.i,ind.j,rinbr).B * fA; 

  double muK = coeff_mu(aK, aL); 
  double muL = 1.0 - muK; 

  double B_sigma = muL * aL - muK * aK;
  double B_sigma_p = bsigma_plus(B_sigma);
  double B_sigma_n = bsigma_minus(B_sigma);

  double A_K = muK * (alpha_osf_(i,j, inbr).A + alpha_osf_(i,j, inbr).B) + B_sigma_p / (f_(i, j) + 1e-15);
  double A_L = muL * (alpha_osf_(ind.i,ind.j, rinbr).A + alpha_osf_(ind.i,ind.j, inbr).B) + B_sigma_n / (f_(ind.i, ind.j) + 1e-15);

  long ii = m.ind2to1(i,j), jj = m.ind2to1(ind.i, ind.j); 
  M_coeffs_.push_back(T(ii, ii, A_K));
  M_coeffs_.push_back(T(ii, jj, -A_L));

}

void Solver::coeff_add_dirbc(int i, int j, int inbr) {  
  Edge edge; 
  m.get_nbr_edg(i, j, inbr, &edge); 

  Ind indA, indB; 
  m.indO(edge.A, &indA); 
  double fA = vertex_f_(indA.i, indA.j); 

  m.indO(edge.B, &indB); 
  double fB = vertex_f_(indB.i,indB.j); 

  double aK = alpha_osf_(i,j,inbr).A * fA + alpha_osf_(i,j,inbr).B * fB; 

  long ii = m.ind2to1(i,j);
  R_(ii) += aK; 
  M_coeffs_.push_back(T(ii, ii, (alpha_osf_(i,j,inbr).A + alpha_osf_(i,j,inbr).B)));
}


void Solver::assemble(){ // obtain M and R 

  // inner cells
  for (std::size_t i=1; i<m.nx()-1; ++i){
    for (std::size_t j=1; j<m.ny()-1; ++j){

      for (size_t inbr = 0; inbr < m.nnbrs(); ++inbr) 
        coeff_add_inner(i, j, inbr);
    }
  }

  // row: i == 0 
  for (std::size_t j=0; j<m.ny(); ++j) {
    // coeff_add_inner(0,j,m.inbr_ip()); // nothing special for alpha=0 for inbr=inbr_im()
    coeff_add_dirbc(0,j,m.inbr_im());  
  }
  for (std::size_t j=1; j<m.ny(); ++j) coeff_add_inner(0,j,m.inbr_jm());
  for (std::size_t j=0; j<m.ny()-1; ++j) coeff_add_inner(0,j,m.inbr_jp());

  // row i == nx-1 alpha = 90 
  for (std::size_t j=0; j<m.ny(); ++j) {
    coeff_add_inner(m.nx()-1,j,m.inbr_im()); 
    // nothing special for alpha=90 for inbr=inbr_ip()
  }
  for (std::size_t j=1; j<m.ny(); ++j) coeff_add_inner(m.nx()-1,j,m.inbr_jm());
  for (std::size_t j=0; j<m.ny()-1; ++j) coeff_add_inner(m.nx()-1,j,m.inbr_jp());

  // col j = 0
  for (std::size_t i=0; i<m.nx(); ++i) {
    coeff_add_inner(i,0,m.inbr_jp()); 
    coeff_add_dirbc(i,0,m.inbr_jm()); 
  }

  for (std::size_t i=1; i<m.nx(); ++i) coeff_add_inner(i,0,m.inbr_im());
  for (std::size_t i=0; i<m.nx()-1; ++i) coeff_add_inner(i,0,m.inbr_ip());

  // col j == ny-1
  for (std::size_t i=0; i<m.nx(); ++i) {
    coeff_add_inner(i,m.ny()-1,m.inbr_jm()); 
    coeff_add_dirbc(i,m.ny()-1,m.inbr_jp()); 
  }

  for (std::size_t i=1; i<m.nx(); ++i) coeff_add_inner(i,m.ny()-1,m.inbr_im());
  for (std::size_t i=0; i<m.nx()-1; ++i) coeff_add_inner(i,m.ny()-1,m.inbr_ip());

  double a, p;
  long ii;
  double Uii;

  for (std::size_t i=0; i<m.nx(); ++i) {
    a = m.x(i);
    for (std::size_t j=0; j<m.ny(); ++j) {
      p = m.y(j);
      ii = m.ind2to1(i,j);
      Uii = G(a, p) * m.area_dt();
      M_coeffs_.push_back(T(ii, ii, Uii));
      R_(ii) += Uii * f_(i,j);
    }
  }

  M_.setFromTriplets(M_coeffs_.begin(), M_coeffs_.end());
}


void Solver::update() {
  R_.setZero();
  M_coeffs_.clear();

  assemble(); 

  solver.analyzePattern(M_);
  solver.factorize(M_);
  f_.reshaped() = solver.solve(R_);
  
  t_ += m.dt(); 
  construct_alpha_osf();
  update_vertex_f();
}

void Solver::update_vertex_f(){
  for (std::size_t i=1; i<m.nx(); ++i)
    for (std::size_t j=1; j<m.ny(); ++j){
      vertex_f_(i,j) = (f_(i-1,j-1) + f_(i-1,j) + f_(i,j-1) + f_(i,j)) / 4.0; 
    }

  double a0, p0;

  // i == 0 and m.nx() boundary
  for (std::size_t j = 1; j<m.ny(); ++j){
    p0 = m.yO() + j*m.dy();
    vertex_f_(0, j) = bcs.alpha0_lc(t(), p0);
    vertex_f_(m.nx(), j) = vertex_f_(m.nx()-1,j);
  }

  // j == 0  and j == m.ny() boundary
  for (std::size_t i = 0; i <= m.nx(); ++i) {
    a0 = m.xO() + i*m.dx(); 
    vertex_f_(i,0) = bcs.pmin(t(), a0); 
    vertex_f_(i,m.ny()) = bcs.pmax(t(), a0); 
  }
}
