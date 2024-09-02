/*
 * File:        D.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "D.h"
#include "common.h"

D::D(const Parameters& paras_in, const Mesh& mesh_in) : paras(paras_in), m(mesh_in) {
    // Initialize matrices Dap, Dpp, and Daa based on mesh size
    Daa_ = Eigen::MatrixXd::Zero(m.nx(), m.ny());
    Dap_ = Eigen::MatrixXd::Zero(m.nx(), m.ny());
    Dpp_ = Eigen::MatrixXd::Zero(m.nx(), m.ny());

    Day_ = Eigen::MatrixXd::Zero(m.nx(), m.ny());
    Dyy_ = Eigen::MatrixXd::Zero(m.nx(), m.ny());

    constructD(paras, 0.0);
}

void D::updateCoefficients(double t) {
    // Implement according to your logic to update Dap, Dpp, Daa with time t
}


// read diffusion coefficients from file
void D::read_d(std::string address, Eigen::MatrixXd* D_rawp){
    Eigen::MatrixXd& D_raw = *D_rawp; 

    std::ifstream fin(address);
    assert(fin.is_open());
    std::string line;

    const double denormalize_factor = gME * gME * gC * gC;
    const double second_to_day = 3600 * 24;

    int nalpha0, nenergy;

    nalpha0 = paras.nalpha0_D();
    nenergy = paras.nE_D();

    for (int i = 0; i < nalpha0; i++){
        for (int j = 0; j < nenergy; j++){
            fin >> D_raw(i, j);
            D_raw(i, j) *= denormalize_factor * second_to_day;
        }
    }
}

void D::locate(double alpha0, double p, Loc* locp){
    int i0, j0;
    double wi, wj;
    double logE = log(p2e(p, gE0)); 

    double pos_alpha0 = (alpha0 - paras.alpha0_min_D()) / paras.dalpha0_D(); 
    double pos_p = (logE - log(paras.Emin_D())) / paras.dlogE_D(); 

    i0 = floor(pos_alpha0); 
    j0 = floor(pos_p); 

    if (i0 >= 0 && i0 < paras.nalpha0_D()) {
      wi = 1 - (pos_alpha0 - i0); 
    } 
    else if (i0 < 0) {
      i0 = 0;
      wi = 1;
    }
    else if (i0 >= paras.nalpha0_D()) {
      i0 = paras.nalpha0_D() - 2; 
      wi = 0.0; 
    }

    if (j0>=0 && j0 < paras.nE_D()) { 
      wj = 1 - (pos_p - j0); 
    }
    else if (j0 < 0) { 
      j0 = 0; 
      wj = 1; 
    }
    else if (j0 >= paras.nE_D()) {
      j0 = paras.nE_D() - 2; 
      wj = 0.0; 
    }

    locp->i0 = i0;
    locp->j0 = j0; 
    locp->wi = wi;
    locp->wj = wj; 
}

double Dinterp(const Eigen::MatrixXd& Draw, const Loc& loc){
   int i0,j0;
   double wi, wj; 
   i0 = loc.i0;
   j0 = loc.j0;
   wi = loc.wi;
   wj = loc.wj; 

   return Draw(i0,j0)*wi*wj + Draw(i0+1,j0)*(1-wi)*wj + Draw(i0+1,j0+1)*(1-wi)*(1-wj) + Draw(i0,j0+1)*wi*(1-wj);
} 

void D::constructD(const Parameters& par, double t){
    double a, p;

    Eigen::MatrixXd Daa_raw(par.nalpha0_D(), par.nE_D());
    Eigen::MatrixXd Dap_raw(par.nalpha0_D(), par.nE_D());
    Eigen::MatrixXd Dpp_raw(par.nalpha0_D(), par.nE_D());

    std::string dfile_base = "D/" + par.dID() + "/" + par.dID() + ".";

    read_d(dfile_base + "Daa", &Daa_raw);
    read_d(dfile_base + "Dap", &Dap_raw);
    read_d(dfile_base + "Dpp", &Dpp_raw);

    Loc loc; 

    for(std::size_t i = 0; i < m.nx(); i++){
        a = m.x(i);
        for(std::size_t j = 0; j < m.ny(); j++){
            p = m.p(j);
            
            locate(a, p, &loc);  
            
            Daa_(i,j) = Dinterp(Daa_raw, loc) / (p*p); 
            Dap_(i,j) = Dinterp(Dap_raw, loc) / p; 
            Dpp_(i,j) = Dinterp(Dpp_raw, loc);

            Day_(i,j) = Dap_(i,j) / p; 
            Dyy_(i,j) = Dpp_(i,j) / (p*p); 
        }
    }
}

