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

struct loc_info{
  double a_dec; // decimal part of coordinate a
  double logE_dec;
  int a_floor; // if a is just on an edge, take opposite a value as a marker
  int logE_floor;
}; 


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

    void read_d(std::string address, Eigen::MatrixXd& D_raw);

    void locate(const Parameters& par, double a, double p, loc_info& loc);
};

#endif /* D_H_ */

