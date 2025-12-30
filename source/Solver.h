/*
 * File:        Solver.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025
 *
 * Copyright (c) Xin Tao
 *
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "common.h"
#include "Mesh.h"
#include "Equation.h"
#include "Parameters.h"

enum class Face { IM, IP, JM, JP };

class Solver {
  public:
    Solver(const Parameters& paras_in, const Mesh& m_in, Equation* eqp);

    void update();
    double t() const { return istep_ * m.dt(); }
    const Xtensor2d& f() const { return f_; }
    double f(const Ind& ind) const { return f_(ind.i, ind.j); }

  private:
    const Parameters& paras;
    const Mesh& m;

    Equation& eq;

    std::size_t istep_; // used to calcualte time = istep_ * dt()

    // M f = R
    Eigen::BiCGSTAB<SpMat> iterSolver;

    SpMat M_;
    std::vector<T> M_coeffs_;

    Xtensor2d f_;
    Eigen::VectorXd R_;
    Eigen::VectorXd ftmp_; // used to store results from the Eigen Solver.

    xt::xtensor<Eigen::Matrix2d, 2> Lambda_;  // Lambda_(i,j) would be the Lambda Matrix at cell (i,j)
    const Eigen::Matrix2d& Lambda(int i, int j) const { return Lambda_(i,j); }
    const Eigen::Matrix2d& Lambda(const Ind& ind) const { return Lambda_(ind.i,ind.j); }
    void update_Lambda();

    //
    // f at vertices to build a lookup table
    //
    Xtensor2d vertex_f_;

    double vertex_f(int i, int j) const { return vertex_f_(i,j); }
    double vertex_f(const Ind& ind) const { return vertex_f_(ind.i, ind.j); }
    void update_vertex_f();



    void assemble();

    // to calculate coefficients alpha_sigma_i and a_sigma_i
    void a_sigma_func(const Ind& ind, int inbr, Vector2* a_sigma_i_p, double* a_sigmap);

    // update coefficients in M for cell Ind, and its neighbor.
    // note that coefficients for both cell Ind and its neighbor are updated.
    void update_coeff_inner_pair(const Ind& ind, int inbr);

    // Here: the inbr neighbor is a Dirichlet boundary.
    void update_coeff_dirbc(const Ind& ind, int inbr);

    void calculate_mu(double a_sigma_K, double a_sigma_L, double* mu_Kp, double* mu_Lp){
      double denom = std::abs(a_sigma_K) + std::abs(a_sigma_L) + 2*gEPS;
      (*mu_Kp) = (std::abs(a_sigma_L)+gEPS)/denom;
      (*mu_Lp) = 1 - (*mu_Kp);
    }

    void init();

    int face_to_inbr(Face face) const {
      switch (face) {
        case Face::IM: return m.inbr_im();
        case Face::IP: return m.inbr_ip();
        case Face::JM: return m.inbr_jm();
        case Face::JP: return m.inbr_jp();
      }
      return m.inbr_im(); // defensive
    }

    Face opposite_face(Face face) const {
      switch (face) {
        case Face::IM: return Face::IP;
        case Face::IP: return Face::IM;
        case Face::JM: return Face::JP;
        case Face::JP: return Face::JM;
      }
      return Face::IM;

    };

    void apply_inner_face_pair(const Ind& cell, Face face) {
      update_coeff_inner_pair(cell, face_to_inbr(face));
    }

    void apply_dirichlet_face(const Ind& cell, Face face) {
      update_coeff_dirbc(cell, face_to_inbr(face));
    }

    // Placeholder for later (currently does nothing)
    void apply_neumann_face(const Ind& cell, Face face) {
      (void)cell;
      (void)face;
      // TODO: implement when BC is modularized
    }

    void apply_boundary_faces_(); 

    void reconstruct_vertex_f_interior_();
    void apply_vertex_bcs_(); // currently eq.apply_bcs(&vertex_f_)

};

#endif /* SOLVER_H */
