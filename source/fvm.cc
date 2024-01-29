/*
 * File:        fvm.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        01/28/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "fvm.h"
#include "Grid.h"

// calculate flux for boundary cells
double flux_bc(const Grid& g, const Matrix& f, int i, int j) {
  return 0;
}

// calculate flux for inner cells
double flux_inner(const Grid& g, const Matrix& f, int i, int j) {

  return 0;
}

void fvm_update(const Grid& g, const Matrix& fn, Matrix* fnp1p) {

  Matrix& fnp1 = *fnp1p; 

  // inner cells
  for (int i=1; i<fnp1.nrows()-1; ++i)
    for (int j=1; j<fnp1.ncols()-1; ++j) {
      fnp1(i,j) = fn(i,j) + g.dt()*flux_inner(g, fn, i, j);
    }


  // boundary cells
  int ibc; 
  for (int j=0; j < fnp1.ncols(); ++j) {
    ibc = 0; 
    fnp1(ibc,j) = fn(ibc,j) + g.dt()*flux_bc(g, fn, ibc, j);

    ibc = fnp1.ncols()-1;
    fnp1(ibc, j) = fn(ibc,j) + g.dt()*flux_bc(g, fn, ibc, j);
  }

  int jbc;
  for (int i=0; i < fnp1.nrows(); ++i) {
    jbc = 0; 
    fnp1(i,jbc) = fn(i,jbc) + g.dt()*flux_bc(g, fn, i, jbc); 

    jbc = fnp1.ncols() - 1; 
    fnp1(i,jbc) = fn(i,jbc) + g.dt()*flux_bc(g, fn, i, jbc); 
  }

}

