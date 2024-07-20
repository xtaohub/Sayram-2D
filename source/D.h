/*
 * File:        D.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef D_H_
#define D_H_

#include "common.h"
#include "Parameters.h"
#include "Mesh.h"

class D {
public:
    D(const Parameters& paras_in, const Mesh& mesh_in);

    double getDap(double t, int i, int j) const;
    double getDpp(double t, int i, int j) const;
    double getDaa(double t, int i, int j) const;

    void constructD(const Parameters& par, double t);

private:
    
    const Parameters& paras; 
    const Mesh& m; 

    Eigen::MatrixXd Daa;
    Eigen::MatrixXd Dap;
    Eigen::MatrixXd Dpp;

    // Update diffusion coefficients with time
    void updateCoefficients(double t);

    void read_d(const Parameters& par, std::string address, Eigen::MatrixXd& D_raw);

    std::vector<double> locate(const Parameters& par, int i, int j);
};

#endif /* D_H_ */

