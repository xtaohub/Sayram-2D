/*
 * File:        Solver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025
 *
 * Copyright (c) Xin Tao
 *
 */

#include "Solver.h"

Solver::Solver(const Parameters& paras_in, const Mesh& m_in, Equation* eqp)
  : paras(paras_in), m(m_in), eq(*eqp){

    std::size_t nx = m.nx();
    std::size_t ny = m.ny();

    M_.resize(nx*ny, nx*ny);

    f_.resize({nx,ny});

    R_.resize(nx*ny);
    ftmp_.resize(nx*ny);

    Lambda_.resize({nx, ny});
    vertex_f_.resize({nx+1,ny+1});

    istep_ = 0;

    init();

  }

void Solver::init(){

  iterSolver.setTolerance(1e-8);

  for (std::size_t i = 0; i < m.nx(); i++){
    for (std::size_t j = 0; j < m.ny(); j++){
      f_(i,j) = eq.init_f({i,j});
    }
  }

  update_Lambda();
  update_vertex_f();
}

void Solver::update_Lambda(){

  for (std::size_t i = 0; i < m.nx(); i++){
    for (std::size_t j = 0; j < m.ny(); j++){
      Lambda_(i,j) << eq.Dxx({i,j}) * eq.G({i,j}), eq.Dxy({i,j}) * eq.G({i,j}),
                      eq.Dxy({i,j}) * eq.G({i,j}), eq.Dyy({i,j}) * eq.G({i,j});
    }
  }
}


void Solver::a_sigma_func(const Ind& ind, int inbr, Vector2* a_sigma_i_p, double* a_sigmap){

  Edge edge;
  m.get_nbr_edge(ind, inbr, &edge);

  Point K = {m.x(ind.i), m.y(ind.j)};

  Point A = edge.v(0);
  Point B = edge.v(1);

  Ind indA, indB;

  Vector2 vkb = B - K;  // as Vector KB (from K point to B)
  Vector2 vka = A - K;

  Vector2 rvkb = {vkb(1), -vkb(0)};  // (x,y) rotated by 90 clockwise, becoming (y, -x)
  Vector2 rvka = {vka(1), -vka(0)};

  double a_sigma_a = edge.length() * (edge.n().transpose() * Lambda(ind) * rvkb)(0) / (vka.transpose() * rvkb)(0);
  double a_sigma_b = edge.length() * (edge.n().transpose() * Lambda(ind) * rvka)(0) / (vkb.transpose() * rvka)(0);

  (*a_sigma_i_p)(0) = a_sigma_a;
  (*a_sigma_i_p)(1) = a_sigma_b;

  m.indO(A, &indA);
  m.indO(B, &indB);

  *a_sigmap = a_sigma_a * vertex_f(indA) + a_sigma_b * vertex_f(indB);
}


void Solver::update_coeff_inner_pair(const Ind& ind, int inbr){ // add coefficient from a inner (no-a-boundary) face

  Ind nbr_ind;

  int rinb;

  double a_sigma_K, a_sigma_L;
  double mu_K, mu_L;

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K, A_L;

  Vector2 a_sigma_i_K;
  Vector2 a_sigma_i_L;

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_i_K, &a_sigma_K);

  // calculate A_L
  m.get_nbr_ind(ind, inbr, &nbr_ind);
  rinb = m.rinbr(inbr);
  a_sigma_func(nbr_ind, rinb, &a_sigma_i_L, &a_sigma_L);

  calculate_mu(a_sigma_K, a_sigma_L, &mu_K, &mu_L);

  B_sigma = mu_L*a_sigma_L - mu_K*a_sigma_K;

  B_sigma_abs = std::abs(B_sigma);
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = mu_K * (a_sigma_i_K[0] + a_sigma_i_K[1]) + B_sigma_plus/(f(ind)+gEPS);
  A_L = mu_L * (a_sigma_i_L[0] + a_sigma_i_L[1]) + B_sigma_minus/(f(nbr_ind)+gEPS);

  std::size_t K = m.flatten_index(ind), L = m.flatten_index(nbr_ind);

  M_coeffs_.push_back(T(K, K, A_K));
  M_coeffs_.push_back(T(K, L, -A_L));
  M_coeffs_.push_back(T(L, L, A_L));
  M_coeffs_.push_back(T(L, K, -A_K));

}

void Solver::update_coeff_dirbc(const Ind& ind, int inbr){ // add coefficient from a Dirichlet face

  double a_sigma_K;

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K;
  Vector2 a_sigma_i_K;

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_i_K, &a_sigma_K);

  B_sigma = -a_sigma_K;
  B_sigma_abs = std::abs(B_sigma);
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = a_sigma_i_K[0] + a_sigma_i_K[1] + B_sigma_plus/(f(ind)+gEPS);

  std::size_t K = m.flatten_index(ind);
  M_coeffs_.push_back(T(K, K, A_K));
  R_(K) += B_sigma_minus;
}


void Solver::assemble(){ // obtain M and R

  std::size_t i, j;

  // i direction
  // i==0: dirichlet boundary condition

  i = 0;
  for (j=0; j<m.ny(); ++j)
      update_coeff_dirbc({i,j}, m.inbr_im());

  for (i=1; i<m.nx(); ++i)
    for (j=0; j<m.ny(); ++j)
          update_coeff_inner_pair({i,j}, m.inbr_im());

// i==m.nx()-1;
// put the handling of the boundary conditions at i==nx-1 here.
// it's neumann type (df/dx=0) in our case, so do nothing.

  // j direction
  // j == 0: handle dirbc
  j = 0;
  for (i=0; i<m.nx(); ++i)
    update_coeff_dirbc({i,j}, m.inbr_jm());

  for (i=0; i<m.nx(); ++i)
    for (j=1; j<m.ny(); ++j)
      update_coeff_inner_pair({i,j}, m.inbr_jm());

  // j == ny-1: handle dirbc
  j = m.ny()-1;
  for (i=0; i<m.nx(); ++i)
    update_coeff_dirbc({i,j}, m.inbr_jp());

  long K;
  double UKK;

  for (std::size_t i=0; i<m.nx(); ++i) {
    for (std::size_t j=0; j<m.ny(); ++j) {
      K = m.flatten_index({i,j});
      UKK = eq.G({i,j}) * m.cell_area_dt();
      M_coeffs_.push_back(T(K, K, UKK));
      R_(K) += UKK * f_(i,j);
    }
  }

  M_.setFromTriplets(M_coeffs_.begin(), M_coeffs_.end());
}


void Solver::update() {
  R_.setZero();
  M_coeffs_.clear();

  assemble();

  iterSolver.compute(M_);
  ftmp_ = iterSolver.solve(R_);

  for (std::size_t i=0; i<m.nx(); ++i){
    for (std::size_t j=0; j<m.ny(); ++j) {
      f_(i,j) = ftmp_(m.flatten_index({i, j}));
    }
  }

  istep_ += 1;
  eq.update(t());
  update_Lambda();
  update_vertex_f();
}

void Solver::update_vertex_f(){

  for (std::size_t i=1; i<m.nx(); ++i){
    for (std::size_t j=1; j<m.ny(); ++j){
      vertex_f_(i,j) = (f_(i-1,j-1) + f_(i-1,j) + f_(i,j-1) + f_(i,j)) / 4.0;
    }
  }

  eq.apply_bcs(&vertex_f_);

}
