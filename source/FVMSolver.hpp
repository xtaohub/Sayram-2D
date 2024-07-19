/*
 * File:        FVMsolver.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef FVM_SOLVER_H
#define FVM_SOLVER_H

// #include "mkl.h"
#include "Mesh.h"
#include "D.h"
#include "BoundaryConditions.h"
#include "Eigen/Core"
#include "Parameters.hpp"
#include "Eigen/Sparse"
#include "Array.h"

struct NTPFA_nodes{ // the two points A,B used in Nonlinear Two Point Approximation
  double A;
  double B;
}; 

// struct Alpha_K{ // The AlphaK matrix for each cell with four faces
//   // NTPFA_nodes ip1;
//   // NTPFA_nodes im1;
//   // NTPFA_nodes jp1;
//   // NTPFA_nodes jm1;
//   double A; 
//   double B; 
// };

class FVMSolver {
  public:
    FVMSolver(const Mesh& m_in, const D& d_in, const BoundaryConditions& bc_in);

    void update();
    const Eigen::MatrixXd& f() const { return f_; }

    void initial();

  private:
    const Mesh& m;
    const D& d; 
    const BoundaryConditions& bc;

    // M f = S
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

    SpMat V_invU_;
    std::vector<T> V_invU_coeffs_;
    SpMat M_;

    Eigen::MatrixXd f_; 
    Eigen::VectorXd S_invU_;
    Eigen::VectorXd R_;

    Eigen::MatrixXd Id_; 

    // Eigen::Matrix<Alpha_K, Eigen::Dynamic, Eigen::Dynamic> alpha_K_;

    Array<NTPFA_nodes, 3> alpha_K_; 

    double hdx_; 
    double hdy_; 

    //
    // use a matrix to store f at vertices to build a lookup table for fA and fB
    // vertex_f is of size (nx+1, ny+1)
    // 
    Array<double, 2> vertex_f_; 

    int ind2to1(int i, int j) const { 
      // map 2d indices to 1, column major
      return j*m.nx()+i; 
    }

    double vertex_f(int imin, int jmin) const { 
      // 
      // obtain flux at vertex by four-point inteprolation.
      // here (i,j) is the smallest indices among the four cells.
      // 
      //          | 
      //   (i,j)  |  (i,j+1)  
      //  --------v-----------
      //  (i+1,j) |  (i+1,j+1)
      //          |
      //
      return (f_(imin,jmin) + f_(imin,jmin+1) + f_(imin+1,jmin) + f_(imin+1,jmin+1)) / 4.0;
    }

    void set_vertex_f(); 

    void assemble();
    void construct_alpha_K();

    void alpha_K(const Eigen::Matrix2d& Lambda_K, const Eigen::Vector2d& K, const Eigen::Vector2d& A, const Eigen::Vector2d& B, double* alpha_KA, double* alpha_KB); 

    // // add coefficients to the Matrix M_ for inner grids (j plus (jp) direction flux)
    // void coeff_add_jp(int i, int j, double a, double p, double u_ipjp, double u_imjp);
    //
    // // add coefficients to the Matrix M_ for inner grids (jm)
    // void coeff_add_jm(int i, int j, double a, double p, double u_imjm, double u_ipjm);
    //
    // // add coefficients to the Matrix M_ for inner grids (ip)
    // void coeff_add_ip(int i, int j, double a, double p, double u_ipjm, double u_ipjp);
    //
    // // add coefficients to the Matrix M_ for inner grids (im)
    // void coeff_add_im(int i, int j, double a, double p, double u_imjp, double u_imjm);

    // add coefficients, but the edge condition (i max boundary condition flux is 0, ignored)
    void coeff_add_jp_edge(int i, int j, double a, double p, double u_ipjp, double u_imjp);
    void coeff_add_jm_edge(int i, int j, double a, double p, double u_imjm, double u_ipjm);
    void coeff_add_im_edge(int i, int j, double a, double p, double u_imjp, double u_imjm);

    void coeff_M(int icell, int jcell);  

    // void coeff_M_single(int i, int j, double a, double p, double u1, double u2); 
    
    // add coefficient to M corresponds to the inbr cell of cell (i,j)
    void coeff_M_add_inner(int i, int j, int inbr);  

    double coeff_a(double alphaA, double fA, double alphaB, double fB) const {
      return alphaA*fA + alphaB*fB; 
    }

    double coeff_mu(double aK, double aL) const {
      if (aK != 0 || aL != 0){
        return abs(aL) / (abs(aK) + abs(aL));
      } else {
        return 0.5;
      }
    }

    double bsigma_inner_cells(double aK, double muK, double aL, double muL) const { 
      return muL*aL-muK*aK;  
    }

    double bsigma_boundary_cells(double aK){
      return -aK; 
    }

    double bsigma_plus(double bsigma){
      return (std::abs(bsigma) + std::abs(bsigma))/2.0;
    }

    double bsigma_minus(double bsigma){
      return (std::abs(bsigma) - std::abs(bsigma))/2.0;
    }

    double init_f(double a, double p){
      return exp(-(p2e(p, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (p * p);
    }

    double boundary_imin(double p){
        return 0.0;
    }

    double boundary_jmin(double a){
      return exp(-(p2e(P_MIN, gE0) - 0.2) / 0.1) * (sin(a) - sin(ALPHA_LC)) / (P_MIN * P_MIN);
    }

    double boundary_jmax(double a){
      return 0.0;
    }

};

#endif /* FVM_SOLVER_H */
