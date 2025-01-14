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

#include <vector>
#include "common.h"
#include "Parameters.h" 
#include "Edge.h"

struct Ind{
  std::size_t i;
  std::size_t j;
}; 

class Mesh {
  public:
    Mesh(const Parameters& p); 

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

    double cell_area() const { return dx() * dy(); }
    double cell_area_dt() const { return cell_area() / dt();}

    void indO(const Point& A, Ind* indp) const { 
      // calculate the i,j coordinate relative to the Origin
      // Note: not the cell index.
      // This function is useful to calculate fA and fB from interpolation
      indp->i = round((A(0) - xO()) / dx());
      indp->j = round((A(1) - yO()) / dy()); 
      
    }

    std::size_t flatten_index(const Ind& ind) const { // map 2d indices to 1
      return ind.j*nx() + ind.i;
    }

    std::size_t nnbrs() const { return 4; } // each cell has 4 nbrs

    // define the neighbor # of 4 adjacent cells
    // im -- (i-1, j, k); jp -- (i, j+1, k)
    // ip -- (i+1, j, k); jm -- (i, j-1, k)
    int inbr_im() const { return 0; }
    int inbr_jp() const { return 1; }
    int inbr_ip() const { return 2; }
    int inbr_jm() const { return 3; }

    int rinbr(int inbr) const { return rinbr_(inbr); }                                    

    void get_nbr_ind(const Ind& ind, int inbr, Ind* nbr_indp) const {
      *nbr_indp = nbr_inds(ind.i,ind.j,inbr); 
    }

    void get_nbr_edge(const Ind& ind, int inbr, Edge* edgep) const {
      *edgep = edges(ind.i,ind.j,inbr); 
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

    Eigen::Matrix<int, 4, 1> rinbr_; 

    xt::xtensor<Ind,3> nbr_inds;
    xt::xtensor<Edge,3> edges;

    void build_connectivity();
};

#endif /* MESH_H */

