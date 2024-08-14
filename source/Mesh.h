/*
 * File:        Mesh.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef MESH_H_
#define MESH_H_

#include "common.h"
#include "Parameters.h" 
#include <vector>
#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>

struct Ind{
  int i; 
  int j;
}; 

typedef Eigen::Vector2d Point;  // each point has two coordinates, Point(0) -- x, Point(1) -- y

struct Edge{
  Point A = {0, 0}; 
  Point B = {0, 0};
}; 

class Mesh {
  public:
    Mesh(const Parameters& p): x_(p.nalpha0()), y_(p.nE()) {

        nx_ = p.nalpha0(); 
        ny_ = p.nE();
        dt_ = p.dt(); 

        xO_ = p.alpha0_lc(); 
        yO_ = p.pmin();

        dx_ = (p.alpha0_max() - p.alpha0_lc()) / p.nalpha0(); 
        dy_ = (p.pmax() - p.pmin()) / p.nE(); 

        x_(0) = xO() + dx()/2.0; 
        y_(0) = yO() + dy()/2.0; 

        for (std::size_t i=1; i<nx(); ++i) x_(i) = x_(0) + i*dx(); 
        for (std::size_t j=1; j<ny(); ++j) y_(j) = y_(0) + j*dy(); 

        nbr_inds.resize({nx(), ny(), 4});
        edges.resize({nx(), ny(), 4});

        build_connectivity(); 

        // The reverse inbr number.
        // For example, if the current cell is K, its 0th neighbor is L.
        // Then, for cell L, K is its 2th neighbor.
        rinbr_(0) = 2; 
        rinbr_(1) = 3; 
        rinbr_(2) = 0; 
        rinbr_(3) = 1; 

      }

    const Eigen::VectorXd& x() const { return x_; }
    const Eigen::VectorXd& y() const { return y_; }

    double x(int i) const { return x_(i); }
    double y(int j) const { return y_(j); }
    double xO() const { return xO_; }
    double yO() const { return yO_; }
    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double dt() const { return dt_; }
    double area_dt() const { return dx_ * dy_ / dt_;}

    int ind2to1(int i, int j) const { // map 2d indices to 1, column major
      return j*nx()+i; 
    }

    void indO(const Point& A, Ind* indp) const { 
      // calculate the i,j coordinate relative to the Origin
      // Note: not the cell index.
      // This function is useful to calculate fA and fB from interpolation
      indp->i = round((A(0) - xO()) / dx());
      indp->j = round((A(1) - yO()) / dy()); 
    }

    std::size_t nnbrs() const { return 4; } // each cell has 4 nbrs
                                    
    // define the neighbor # of four adjacent cells
    // im -- (i-1, j); jp -- (i,j+1)
    // ip -- (i+1, j); jm -- (i,j-1)
    int inbr_im() const { return 0; }
    int inbr_jp() const { return 1; }
    int inbr_ip() const { return 2; }
    int inbr_jm() const { return 3; }

    int rinbr(int inbr) const { return rinbr_(inbr); }                                    

    void get_nbr_ind(int i, int j, int inbr, Ind* nbr_indp) const {
      *nbr_indp = nbr_inds(i,j,inbr); 
    }

    void get_nbr_edg(int i, int j, int inbr, Edge* edgep) const {
      *edgep = edges(i,j,inbr); 
    }

  private:
    std::size_t nx_;
    std::size_t ny_;
    double dx_;
    double dy_;
    double dt_; 

    // coordinate origin: corresponds to i-0.5, j-0.5
    double xO_; 
    double yO_; 

    Eigen::VectorXd x_; 
    Eigen::VectorXd y_; 

    Eigen::Vector4i rinbr_; 

    xt::xtensor<Ind,3> nbr_inds;
    xt::xtensor<Edge,3> edges;


    void build_connectivity();
};

#endif /* MESH_H */

