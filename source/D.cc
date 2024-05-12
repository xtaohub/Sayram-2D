/*
 * File:        D.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "D.h"

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

