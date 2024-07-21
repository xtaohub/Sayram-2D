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

// struct loc_info{
//   double a_dec; // decimal part of coordinate a
//   double logE_dec;
//   int a_floor; // if a is just on an edge, take opposite a value as a marker
//   int logE_floor;
// }; 

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
    void locate(double alpha, double p, Loc* locp);
    void read_d(std::string address, Eigen::MatrixXd* D_rawp);
};

#endif /* D_H_ */

