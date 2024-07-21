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

    int nalpha, nenergy;

    nalpha = paras.nalpha_D();
    nenergy = paras.nE_D();

    for (int i = 0; i < nalpha; i++){
        for (int j = 0; j < nenergy; j++){
            fin >> D_raw(i, j);
            D_raw(i, j) *= denormalize_factor * second_to_day;
        }
    }
}

//
// void D::locate(const Parameters& par, double a, double p, loc_info& loc){
//     std::vector<double> locinfo = {0.0, 0.0, 0.0, 0.0};
//     double a_ = (a - par.alpha_min_D()) / gPI * 180.0;
//     double a_f = floor(a_);
//     double logE = (log(p2e(p, gE0)) - log(par.Emin_D())) / par.dlogE_D();
//     double logE_f = floor(logE);
//     loc.a_dec = a_ - a_f;
//     loc.logE_dec = logE - logE_f;
//     // edge judege:
//     if (abs(loc.a_dec) < 1e-6) {
//         loc.a_floor = -round(a_f);
//     } else {
//         loc.a_floor = round(a_f);
//     }
//     if (abs(loc.logE_dec) < 1e-6) {
//         loc.logE_floor = -round(logE_f);
//     } else {
//         loc.logE_floor = round(logE_f);
//     }
// }

void D::locate(double alpha, double p, Loc* locp){
    int i0, j0;
    double wi, wj;
    double logE = log(p2e(p, gE0)); 

    double pos_alpha = (alpha - paras.alpha_min_D()) / paras.dalpha_D(); 
    double pos_p = (logE - log(paras.Emin_D())) / paras.dlogE_D(); 

    i0 = floor(pos_alpha); 
    j0 = floor(pos_p); 

    if (i0 >= 0 && i0 < paras.nalpha_D()) {
      wi = 1 - (pos_alpha - i0); 
    } 
    else if (i0 < 0) {
      i0 = 0;
      wi = 1;
    }
    else if (i0 >= paras.nalpha_D()) {
      i0 = paras.nalpha_D() - 2; 
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

    Eigen::MatrixXd Daa_raw(par.nalpha_D(), par.nE_D());
    Eigen::MatrixXd Dap_raw(par.nalpha_D(), par.nE_D());
    Eigen::MatrixXd Dpp_raw(par.nalpha_D(), par.nE_D());

    std::string dfile_base = "D/" + par.dID() + "/" + par.dID() + ".";

    read_d(dfile_base + "Daa", &Daa_raw);
    read_d(dfile_base + "Dap", &Dap_raw);
    read_d(dfile_base + "Dpp", &Dpp_raw);

    Loc loc; 

    for(int i = 0; i < m.nx(); i++){
        a = m.x(i);
        for(int j = 0; j < m.ny(); j++){
            p = m.y(j);
            
            locate(a, p, &loc);             
            
            Daa_(i,j) = Dinterp(Daa_raw, loc) / (p*p); 
            Dap_(i,j) = Dinterp(Dap_raw, loc) / p; 
            Dpp_(i,j) = Dinterp(Dpp_raw, loc);


            // locate(par, a, p, loc);

            // if(loc.a_floor < 0 and loc.logE_floor < 0){
            //     Daa_(i, j) = Daa_raw(-loc.a_floor, -loc.logE_floor) / (p * p);
            //     Dap_(i, j) = Dap_raw(-loc.a_floor, -loc.logE_floor) / p;
            //     Dpp_(i, j) = Dpp_raw(-loc.a_floor, -loc.logE_floor);
            // }
            // else if(loc.a_floor < 0) {
            //     double w1 = 1.0 / loc.logE_dec;
            //     double w2 = 1.0 / (1.0 - loc.logE_dec);
            //     Daa_(i, j) = (w1 * Daa_raw(-loc.a_floor, loc.logE_floor) + w2 * Daa_raw(-loc.a_floor, loc.logE_floor + 1)) / (w1 + w2) / (p * p);
            //     Dap_(i, j) = (w1 * Dap_raw(-loc.a_floor, loc.logE_floor) + w2 * Dap_raw(-loc.a_floor, loc.logE_floor + 1)) / (w1 + w2) / p;
            //     Dpp_(i, j) = (w1 * Dpp_raw(-loc.a_floor, loc.logE_floor) + w2 * Dpp_raw(-loc.a_floor, loc.logE_floor + 1)) / (w1 + w2);
            // }
            // else if(loc.logE_floor < 0) {
            //     double w1 = 1.0 / loc.a_dec;
            //     double w2 = 1.0 / (1.0 - loc.a_dec);
            //     Daa_(i, j) = (w1 * Daa_raw(loc.a_floor, -loc.logE_floor) + w2 * Daa_raw(loc.a_floor + 1, -loc.logE_floor)) / (w1 + w2) / (p * p);
            //     Dap_(i, j) = (w1 * Dap_raw(loc.a_floor, -loc.logE_floor) + w2 * Dap_raw(loc.a_floor + 1, -loc.logE_floor)) / (w1 + w2) / p;
            //     Dpp_(i, j) = (w1 * Dpp_raw(loc.a_floor, -loc.logE_floor) + w2 * Dpp_raw(loc.a_floor + 1, -loc.logE_floor)) / (w1 + w2);
            // }
            // else {
            //     double w1 = 1.0 / sqrt(loc.a_dec * loc.a_dec + loc.logE_dec * loc.logE_dec);
            //     double w2 = 1.0 / sqrt(loc.a_dec * loc.a_dec + (1.0 - loc.logE_dec) * (1.0 - loc.logE_dec));
            //     double w3 = 1.0 / sqrt((1.0 - loc.a_dec) * (1.0 - loc.a_dec) + (1.0 - loc.logE_dec) * (1.0 - loc.logE_dec));
            //     double w4 = 1.0 / sqrt((1.0 - loc.a_dec) * (1.0 - loc.a_dec) + loc.logE_dec * loc.logE_dec);
            //     Daa_(i, j) = (w1 * Daa_raw(loc.a_floor, loc.logE_floor) + w2 * Daa_raw(loc.a_floor, loc.logE_floor + 1) + w3 * Daa_raw(loc.a_floor + 1, loc.logE_floor + 1) + w4 * Daa_raw(loc.a_floor + 1, loc.logE_floor)) / (w1 + w2 + w3 + w4) / (p * p);
            //     Dap_(i, j) = (w1 * Dap_raw(loc.a_floor, loc.logE_floor) + w2 * Dap_raw(loc.a_floor, loc.logE_floor + 1) + w3 * Dap_raw(loc.a_floor + 1, loc.logE_floor + 1) + w4 * Dap_raw(loc.a_floor + 1, loc.logE_floor)) / (w1 + w2 + w3 + w4) / p;
            //     Dpp_(i, j) = (w1 * Dpp_raw(loc.a_floor, loc.logE_floor) + w2 * Dpp_raw(loc.a_floor, loc.logE_floor + 1) + w3 * Dpp_raw(loc.a_floor + 1, loc.logE_floor + 1) + w4 * Dpp_raw(loc.a_floor + 1, loc.logE_floor)) / (w1 + w2 + w3 + w4);
            // }
        }
    }
}

