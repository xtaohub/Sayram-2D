/*
 * File:        D.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
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
#include "Parameters.hpp"
#include "common.hpp"

D::D(const Grid& grid) : grid(grid) {
    // Initialize matrices Dap, Dpp, and Daa based on grid size
    Daa = Eigen::MatrixXd::Zero(grid.nx(), grid.ny());
    Dap = Eigen::MatrixXd::Zero(grid.nx(), grid.ny());
    Dpp = Eigen::MatrixXd::Zero(grid.nx(), grid.ny());
}

D::D(const Grid& grid, const Eigen::MatrixXd& Daa_in, const Eigen::MatrixXd& Dap_in, const Eigen::MatrixXd& Dpp_in): grid(grid) {

    Daa = Eigen::MatrixXd::Zero(grid.nx(), grid.ny());
    Dap = Eigen::MatrixXd::Zero(grid.nx(), grid.ny());
    Dpp = Eigen::MatrixXd::Zero(grid.nx(), grid.ny());

    Daa = Daa_in;
    Dap = Dap_in;
    Dpp = Dpp_in; 

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
void read_d(std::string address, Eigen::MatrixXd& D_raw){

    std::ifstream fin(address);
    assert(fin.is_open());
    std::string line;

    for (int i = 0; i < 91; i++){
        getline(fin, line);
        if (i > 0){
            for (int j = 0; j < 49; j++){
                try{
                    double temp = std::stod(line.substr((5 + 15) * j + 4, 16));
                    D_raw(j, i - 1) = temp;
                } catch (const std::exception& e){
                    std::cout<< "File "<< address << " has data error." << std::endl;
                }
            }
        }
    }
}


std::vector<double> locate(int j, int i){
    std::vector<double> locinfo = {0.0, 0.0, 0.0, 0.0};
    double a = (ALPHA_LC + hdx + dx * i) / gPI * 180.0 - 1;
    double af = floor(a);
    double p = P_MIN + hdy + dy * j;
    double logE = (log(p2e(p, gE0)) - log(E_RANGE[0])) / dlogE;
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


void D::constructD(double t){
    Eigen::MatrixXd Daa_raw(49, 90);
    Eigen::MatrixXd Dap_raw(49, 90);
    Eigen::MatrixXd Dpp_raw(49, 90);

    std::string address1 = "D/AlbertYoung_chorus/AlbertYoung_chorus.Daa";
    std::string address2 = "D/AlbertYoung_chorus/AlbertYoung_chorus.Dap";
    std::string address3 = "D/AlbertYoung_chorus/AlbertYoung_chorus.Dpp";

    const double denormalize_factor = gME * gME * gC * gC;

    read_d(address1, Daa_raw);
    read_d(address2, Dap_raw);
    read_d(address3, Dpp_raw);

    for(int j = 0; j < ny; j++){
        for(int i = 0; i < nx; i++){
            double p = P_MIN + hdy + dy * j;
            std::vector<double> loc = locate(j, i);
            int locx = int(loc[0]);
            int locy = int(loc[1]);
            double da = loc[2];
            double dp = loc[3];
            if(loc[0] < 0 and loc[1] < 0){
                Daa(j, i) = Daa_raw(-locy, -locx) * denormalize_factor / (p * p);
                Dap(j, i) = Dap_raw(-locy, -locx) * denormalize_factor / p;
                Dpp(j, i) = Dpp_raw(-locy, -locx) * denormalize_factor;
            }
            else if(loc[1] < 0) {
                double w1 = 1.0 / da;
                double w2 = 1.0 / (1.0 - da);
                Daa(j, i) = (w1 * Daa_raw(-locy, locx) + w2 * Daa_raw(-locy, locx + 1)) / (w1 + w2) * denormalize_factor / (p * p);
                Dap(j, i) = (w1 * Dap_raw(-locy, locx) + w2 * Dap_raw(-locy, locx + 1)) / (w1 + w2) * denormalize_factor / p;
                Dpp(j, i) = (w1 * Dpp_raw(-locy, locx) + w2 * Dpp_raw(-locy, locx + 1)) / (w1 + w2) * denormalize_factor;
            }
            else if(loc[0] < 0) {
                double w1 = 1.0 / dp;
                double w2 = 1.0 / (1.0 - dp);
                Daa(j, i) = (w1 * Daa_raw(locy, -locx) + w2 * Daa_raw(locy + 1, -locx)) / (w1 + w2) * denormalize_factor / (p * p);
                Dap(j, i) = (w1 * Dap_raw(locy, -locx) + w2 * Dap_raw(locy + 1, -locx)) / (w1 + w2) * denormalize_factor / p;
                Dpp(j, i) = (w1 * Dpp_raw(locy, -locx) + w2 * Dpp_raw(locy + 1, -locx)) / (w1 + w2) * denormalize_factor;
            }
            else {
                double w1 = 1.0 / sqrt(da * da + dp * dp);
                double w2 = 1.0 / sqrt((1.0 - da) * (1.0 - da) + dp * dp);
                double w3 = 1.0 / sqrt((1.0 - da) * (1.0 - da) + (1.0 - dp) * (1.0 - dp));
                double w4 = 1.0 / sqrt(da * da + (1.0 - dp) * (1.0 - dp));
                Daa(j, i) = (w1 * Daa_raw(locy, locx) + w2 * Daa_raw(locy, locx + 1) + w3 * Daa_raw(locy + 1, locx + 1) + w4 * Daa_raw(locy + 1, locx)) / (w1 + w2 + w3 + w4) * denormalize_factor / (p * p);
                Dap(j, i) = (w1 * Dap_raw(locy, locx) + w2 * Dap_raw(locy, locx + 1) + w3 * Dap_raw(locy + 1, locx + 1) + w4 * Dap_raw(locy + 1, locx)) / (w1 + w2 + w3 + w4) * denormalize_factor / p;
                Dpp(j, i) = (w1 * Dpp_raw(locy, locx) + w2 * Dpp_raw(locy, locx + 1) + w3 * Dpp_raw(locy + 1, locx + 1) + w4 * Dpp_raw(locy + 1, locx)) / (w1 + w2 + w3 + w4) * denormalize_factor;
            }
        }
    }
}

