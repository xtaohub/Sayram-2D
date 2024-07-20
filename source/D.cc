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
void D::read_d(const Parameters& par, std::string address, Eigen::MatrixXd& D_raw){

    std::ifstream fin(address);
    assert(fin.is_open());
    std::string line;

    for (int i = 0; i < par.nalpha_D() + 1; i++){
        getline(fin, line);
        if (i > 0){
            for (int j = 0; j < par.nE_D(); j++){
                try{
                    double temp = std::stod(line.substr((5 + 15) * j + 4, 16));
                    D_raw(i - 1, j) = temp;
                } catch (const std::exception& e){
                    std::cout<< "File "<< address << " has data error." << std::endl;
                }
            }
        }
    }
}


std::vector<double> D::locate(const Parameters& par, int i, int j){
    std::vector<double> locinfo = {0.0, 0.0, 0.0, 0.0};
    double a = (m.x(i)) / gPI * 180.0 - 1;
    double af = floor(a);
    double p = m.y(j);
    double logE = (log(p2e(p, gE0)) - log(par.Emin_D())) / par.dlogE_D();
    double logEf = floor(logE);
    locinfo[2] = a - af;
    locinfo[3] = logE - logEf;
    // edge judege:
    if (abs(locinfo[2]) < 1e-6) {
        locinfo[0] = -af;
    } else {
        locinfo[0] = af;
    }
    if (abs(locinfo[3]) < 1e-6) {
        locinfo[1] = -logEf;
    } else {
        locinfo[1] = logEf;
    }
    return locinfo;
}


void D::constructD(const Parameters& par, double t){
    Eigen::MatrixXd Daa_raw(par.nalpha_D(), par.nE_D());
    Eigen::MatrixXd Dap_raw(par.nalpha_D(), par.nE_D());
    Eigen::MatrixXd Dpp_raw(par.nalpha_D(), par.nE_D());

    std::string address1 = "D/AlbertYoung_chorus/AlbertYoung_chorus.Daa";
    std::string address2 = "D/AlbertYoung_chorus/AlbertYoung_chorus.Dap";
    std::string address3 = "D/AlbertYoung_chorus/AlbertYoung_chorus.Dpp";

    const double denormalize_factor = gME * gME * gC * gC;

    read_d(par, address1, Daa_raw);
    read_d(par, address2, Dap_raw);
    read_d(par, address3, Dpp_raw);

    for(int i = 0; i < m.nx(); i++){
        for(int j = 0; j < m.ny(); j++){
            double p = m.y(j);
            std::vector<double> loc = locate(par, i, j);
            int locx = int(loc[0]);
            int locy = int(loc[1]);
            double da = loc[2];
            double dp = loc[3];
            if(loc[0] < 0 and loc[1] < 0){
                Daa(i, j) = Daa_raw(-locx, -locy) * denormalize_factor / (p * p);
                Dap(i, j) = Dap_raw(-locx, -locy) * denormalize_factor / p;
                Dpp(i, j) = Dpp_raw(-locx, -locy) * denormalize_factor;
            }
            else if(loc[0] < 0) {
                double w1 = 1.0 / dp;
                double w2 = 1.0 / (1.0 - dp);
                Daa(i, j) = (w1 * Daa_raw(-locx, locy) + w2 * Daa_raw(-locx, locy + 1)) / (w1 + w2) * denormalize_factor / (p * p);
                Dap(i, j) = (w1 * Dap_raw(-locx, locy) + w2 * Dap_raw(-locx, locy + 1)) / (w1 + w2) * denormalize_factor / p;
                Dpp(i, j) = (w1 * Dpp_raw(-locx, locy) + w2 * Dpp_raw(-locx, locy + 1)) / (w1 + w2) * denormalize_factor;
            }
            else if(loc[1] < 0) {
                double w1 = 1.0 / da;
                double w2 = 1.0 / (1.0 - da);
                Daa(i, j) = (w1 * Daa_raw(locx, -locy) + w2 * Daa_raw(locx + 1, -locy)) / (w1 + w2) * denormalize_factor / (p * p);
                Dap(i, j) = (w1 * Dap_raw(locx, -locy) + w2 * Dap_raw(locx + 1, -locy)) / (w1 + w2) * denormalize_factor / p;
                Dpp(i, j) = (w1 * Dpp_raw(locx, -locy) + w2 * Dpp_raw(locx + 1, -locy)) / (w1 + w2) * denormalize_factor;
            }
            else {
                double w1 = 1.0 / sqrt(da * da + dp * dp);
                double w2 = 1.0 / sqrt(da * da + (1.0 - dp) * (1.0 - dp));
                double w3 = 1.0 / sqrt((1.0 - da) * (1.0 - da) + (1.0 - dp) * (1.0 - dp));
                double w4 = 1.0 / sqrt((1.0 - da) * (1.0 - da) + dp * dp);
                Daa(i, j) = (w1 * Daa_raw(locx, locy) + w2 * Daa_raw(locx, locy + 1) + w3 * Daa_raw(locx + 1, locy + 1) + w4 * Daa_raw(locx + 1, locy)) / (w1 + w2 + w3 + w4) * denormalize_factor / (p * p);
                Dap(i, j) = (w1 * Dap_raw(locx, locy) + w2 * Dap_raw(locx, locy + 1) + w3 * Dap_raw(locx + 1, locy + 1) + w4 * Dap_raw(locx + 1, locy)) / (w1 + w2 + w3 + w4) * denormalize_factor / p;
                Dpp(i, j) = (w1 * Dpp_raw(locx, locy) + w2 * Dpp_raw(locx, locy + 1) + w3 * Dpp_raw(locx + 1, locy + 1) + w4 * Dpp_raw(locx + 1, locy)) / (w1 + w2 + w3 + w4) * denormalize_factor;
            }
        }
    }
}

