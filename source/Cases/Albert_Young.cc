/*
 * File:        Albert_Young.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Albert_Young.h"


Albert_Young::Albert_Young(const Parameters& paras_in, const Mesh& mesh_in):Equation(mesh_in),paras(paras_in), m(mesh_in), io(paras_in){ // initialize G, D, and etc.
  init();
  constructD(paras, io);
}


double Albert_Young::init_f(const Ind& ind) const { 
  double x, y;
  x = m.x(ind.i);
  y = m.y(ind.j);

  return calculate_init_f(x, y);
}


void Albert_Young::init(){

  double x, y;

  for (std::size_t i = 0; i < m.nx(); i++){
    x = m.x(i);
    for (std::size_t j = 0; j < m.ny(); j++){
      y = m.y(j);
      G_(i,j) = calculate_G(x, y);
    }
  }
}


void Albert_Young::update(double t) {}


void Albert_Young::apply_bcs(Xtensor2d* vertex_fp) const { 
  Xtensor2d& vertex_f = *vertex_fp;  

  double x, y;

  // i == 0 and m.nx() boundary condition case
  for (std::size_t j = 0; j<m.ny()+1; ++j){
    y = m.yO() + j*m.dy();
    vertex_f(0, j) = xmin(y);
    vertex_f(m.nx(), j) = vertex_f(m.nx()-1, j);
  }

  // j == 0  and j == m.ny() boundary
  for (std::size_t i = 0; i<m.nx()+1; ++i){
    x = m.xO() + i*m.dx(); 
    vertex_f(i, 0) = ymin(x); 
    vertex_f(i, m.ny()) = ymax(x);  
  }
}


void Albert_Young::locate(double alpha0, double logE, Loc* locp){
    std::size_t i0, j0;
    double wi, wj;

    double pos_x = (alpha0 - io.xmin_D()) / (io.xmax_D() - io.xmin_D()) * (io.nx_D() - 1);
    double pos_y = (logE - log(io.ymin_D())) / (log(io.ymax_D()) - log(io.ymin_D())) * (io.ny_D() - 1);

    i0 = floor(pos_x); 
    j0 = floor(pos_y); 

    calWeight(i0, wi, io.nx_D()-1, pos_x);
    calWeight(j0, wj, io.ny_D()-1, pos_y);

    locp->i0 = i0;
    locp->j0 = j0; 
    locp->wi = wi;
    locp->wj = wj; 
}

void Albert_Young::constructD(const Parameters& par, const Albert_Young_IO& io){

  const double denormalize_factor = gME * gME * gC * gC;
  const double second_to_day = 3600 * 24;
  double alpha0, logE, p;

  Loc loc; 

  for(std::size_t i = 0; i < m.nx(); i++){
    alpha0 = m.x(i);
    for(std::size_t j = 0; j < m.ny(); j++){
        logE = m.y(j);
        p = e2p(exp(logE), gE0);

        locate(alpha0, logE, &loc);  

        Dxx_(i,j) = interp2D(io.Dxx_raw, loc) * denormalize_factor * second_to_day / (p*p); 
        Dxy_(i,j) = interp2D(io.Dxy_raw, loc) * denormalize_factor * second_to_day * dlogE_dp(logE, gE0) / p;
        Dyy_(i,j) = interp2D(io.Dyy_raw, loc) * denormalize_factor * second_to_day * pow(dlogE_dp(logE, gE0),2); 
      }
  }
}







