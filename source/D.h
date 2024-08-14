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

struct Loc{
  int i0;
  int j0; 
  double wi;
  double wj; 
}; 

class D {
public:
    D(const Parameters& paras_in, const Mesh& mesh_in);

    double Daa(double t, int i, int j) const { return Daa_(i,j); }
    double Dap(double t, int i, int j) const { return Dap_(i,j); }
    double Dpp(double t, int i, int j) const { return Dpp_(i,j); }

    void constructD(const Parameters& par, double t);

private:
    
    const Parameters& paras; 
    const Mesh& m; 

    Eigen::MatrixXd Daa_;
    Eigen::MatrixXd Dap_;
    Eigen::MatrixXd Dpp_;

    // Update diffusion coefficients with time
    void updateCoefficients(double t);
    void locate(double alpha0, double p, Loc* locp);
    void read_d(std::string address, Eigen::MatrixXd* D_rawp);
};

#endif /* D_H_ */

