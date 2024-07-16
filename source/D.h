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

#include "Eigen/Core"
#include "Grid.h"

class D {
public:
    D(const Grid& grid);
    D(const Grid& grid, const Eigen::MatrixXd& Daa_in, const Eigen::MatrixXd& Dap_in, const Eigen::MatrixXd& Dpp_in);

    double getDap(double t, int i, int j) const;
    double getDpp(double t, int i, int j) const;
    double getDaa(double t, int i, int j) const;

    // Update diffusion coefficients with time
    void updateCoefficients(double t);

    void constructD(double t);

private:

    Eigen::MatrixXd Daa;
    Eigen::MatrixXd Dap;
    Eigen::MatrixXd Dpp;

    const Grid& grid;
};

#endif /* D_H_ */

