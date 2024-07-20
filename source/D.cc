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
    Daa = Eigen::MatrixXd::Zero(m.nx(), m.ny());
    Dap = Eigen::MatrixXd::Zero(m.nx(), m.ny());
    Dpp = Eigen::MatrixXd::Zero(m.nx(), m.ny());

    constructD(paras, 0.0);
}

double D::getDap(double t, int i, int j) const {
    // Implement according to your logic
    // Example: return Dap(i, j);
    return Dap(i, j);
}

double D::getDpp(double t, int i, int j) const {
    // Implement according to your logic
    // Example: return Dpp(i, j);
    return Dpp(i, j);
}

double D::getDaa(double t, int i, int j) const {
    // Implement according to your logic
    // Example: return Daa(i, j);
    return Daa(i, j);
}

void D::updateCoefficients(double t) {
    // Implement according to your logic to update Dap, Dpp, Daa with time t
}


// read diffusion coefficients from file
void D::read_d(std::string address, Eigen::MatrixXd& D_raw){

    std::ifstream fin(address);
    assert(fin.is_open());
    std::string line;

    const double denormalize_factor = gME * gME * gC * gC;
    const double second_to_day = 3600 * 24;

    int nalpha, nenergy;
    fin >> nalpha >> nenergy;

    for (int i = 0; i < nalpha; i++){
        for (int j = 0; j < nenergy; j++){
            fin >> D_raw(i, j);
            D_raw(i, j) *= denormalize_factor * second_to_day;
        }
    }
}


void D::locate(const Parameters& par, double a, double p, loc_info& loc){
    std::vector<double> locinfo = {0.0, 0.0, 0.0, 0.0};
    double a_ = (a - par.alpha_min_D()) / gPI * 180.0;
    double a_f = floor(a_);
    double logE = (log(p2e(p, gE0)) - log(par.Emin_D())) / par.dlogE_D();
    double logE_f = floor(logE);
    loc.a_dec = a_ - a_f;
    loc.logE_dec = logE - logE_f;
    // edge judege:
    if (abs(loc.a_dec) < 1e-6) {
        loc.a_floor = -round(a_f);
    } else {
        loc.a_floor = round(a_f);
    }
    if (abs(loc.logE_dec) < 1e-6) {
        loc.logE_floor = -round(logE_f);
    } else {
        loc.logE_floor = round(logE_f);
    }
}


void D::constructD(const Parameters& par, double t){
    double a, p;
    Eigen::MatrixXd Daa_raw(par.nalpha_D(), par.nE_D());
    Eigen::MatrixXd Dap_raw(par.nalpha_D(), par.nE_D());
    Eigen::MatrixXd Dpp_raw(par.nalpha_D(), par.nE_D());

    std::string address1 = "D/" + par.dID() + "/" + par.dID() + ".Daa";
    std::string address2 = "D/" + par.dID() + "/" + par.dID() + ".Dap";
    std::string address3 = "D/" + par.dID() + "/" + par.dID() + ".Dpp";

    read_d(address1, Daa_raw);
    read_d(address2, Dap_raw);
    read_d(address3, Dpp_raw);

    loc_info loc;

    for(int i = 0; i < m.nx(); i++){
        a = m.x(i);
        for(int j = 0; j < m.ny(); j++){
            p = m.y(j);
            locate(par, a, p, loc);
            if(loc.a_floor < 0 and loc.logE_floor < 0){
                Daa(i, j) = Daa_raw(-loc.a_floor, -loc.logE_floor) / (p * p);
                Dap(i, j) = Dap_raw(-loc.a_floor, -loc.logE_floor) / p;
                Dpp(i, j) = Dpp_raw(-loc.a_floor, -loc.logE_floor);
            }
            else if(loc.a_floor < 0) {
                double w1 = 1.0 / loc.logE_dec;
                double w2 = 1.0 / (1.0 - loc.logE_dec);
                Daa(i, j) = (w1 * Daa_raw(-loc.a_floor, loc.logE_floor) + w2 * Daa_raw(-loc.a_floor, loc.logE_floor + 1)) / (w1 + w2) / (p * p);
                Dap(i, j) = (w1 * Dap_raw(-loc.a_floor, loc.logE_floor) + w2 * Dap_raw(-loc.a_floor, loc.logE_floor + 1)) / (w1 + w2) / p;
                Dpp(i, j) = (w1 * Dpp_raw(-loc.a_floor, loc.logE_floor) + w2 * Dpp_raw(-loc.a_floor, loc.logE_floor + 1)) / (w1 + w2);
            }
            else if(loc.logE_floor < 0) {
                double w1 = 1.0 / loc.a_dec;
                double w2 = 1.0 / (1.0 - loc.a_dec);
                Daa(i, j) = (w1 * Daa_raw(loc.a_floor, -loc.logE_floor) + w2 * Daa_raw(loc.a_floor + 1, -loc.logE_floor)) / (w1 + w2) / (p * p);
                Dap(i, j) = (w1 * Dap_raw(loc.a_floor, -loc.logE_floor) + w2 * Dap_raw(loc.a_floor + 1, -loc.logE_floor)) / (w1 + w2) / p;
                Dpp(i, j) = (w1 * Dpp_raw(loc.a_floor, -loc.logE_floor) + w2 * Dpp_raw(loc.a_floor + 1, -loc.logE_floor)) / (w1 + w2);
            }
            else {
                double w1 = 1.0 / sqrt(loc.a_dec * loc.a_dec + loc.logE_dec * loc.logE_dec);
                double w2 = 1.0 / sqrt(loc.a_dec * loc.a_dec + (1.0 - loc.logE_dec) * (1.0 - loc.logE_dec));
                double w3 = 1.0 / sqrt((1.0 - loc.a_dec) * (1.0 - loc.a_dec) + (1.0 - loc.logE_dec) * (1.0 - loc.logE_dec));
                double w4 = 1.0 / sqrt((1.0 - loc.a_dec) * (1.0 - loc.a_dec) + loc.logE_dec * loc.logE_dec);
                Daa(i, j) = (w1 * Daa_raw(loc.a_floor, loc.logE_floor) + w2 * Daa_raw(loc.a_floor, loc.logE_floor + 1) + w3 * Daa_raw(loc.a_floor + 1, loc.logE_floor + 1) + w4 * Daa_raw(loc.a_floor + 1, loc.logE_floor)) / (w1 + w2 + w3 + w4) / (p * p);
                Dap(i, j) = (w1 * Dap_raw(loc.a_floor, loc.logE_floor) + w2 * Dap_raw(loc.a_floor, loc.logE_floor + 1) + w3 * Dap_raw(loc.a_floor + 1, loc.logE_floor + 1) + w4 * Dap_raw(loc.a_floor + 1, loc.logE_floor)) / (w1 + w2 + w3 + w4) / p;
                Dpp(i, j) = (w1 * Dpp_raw(loc.a_floor, loc.logE_floor) + w2 * Dpp_raw(loc.a_floor, loc.logE_floor + 1) + w3 * Dpp_raw(loc.a_floor + 1, loc.logE_floor + 1) + w4 * Dpp_raw(loc.a_floor + 1, loc.logE_floor)) / (w1 + w2 + w3 + w4);
            }
        }
    }
}

