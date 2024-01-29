/*
 * File:        Grid.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        01/28/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef GRID_H_
#define GRID_H_

class Grid {
  public:
    Grid() { 
      dt_ = 1.0;  
      du_ = 1.0; 
      dv_ = 1.0; 
  
      ncu_ = 20;
      ncv_ = 20; 

      G_.resize(ncu(), ncv());
      G_ = 1.0; 
    }; 

    double dt() const { return dt_; }
    double du() const { return du_; }
    double dv() const { return dv_; }

    int ncu() const { return ncu_; }
    int ncv() const { return ncv_; }

    double G(int i, int j) const { return G_(i,j); } // Jacobian
  
  private:
    double dt_; 
    double du_; 
    double dv_; 
    Matrix G_; 
}; 

#endif /* GRID_H */

