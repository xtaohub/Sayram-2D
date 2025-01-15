/*
 * File:        Mesh.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        07/19/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Mesh.h"

Mesh::Mesh(const Parameters& p): x_(p.nalpha0()), y_(p.nE()) {

  nx_ = p.nalpha0(); 
  ny_ = p.nE(); 

  dt_ = p.dt(); 

  xO_ = p.alpha0_min(); 
  yO_ = p.logEmin();

  dx_ = (p.alpha0_max() - p.alpha0_min()) / p.nalpha0(); 
  dy_ = (p.logEmax() - p.logEmin()) / p.nE(); 

  x_(0) = xO() + dx()/2.0; 
  y_(0) = yO() + dy()/2.0; 

  for (std::size_t i=1; i<nx(); ++i) x_(i) = x_(0) + i*dx(); 
  for (std::size_t j=1; j<ny(); ++j) y_(j) = y_(0) + j*dy(); 

  nbr_inds.resize({nx(), ny(), nnbrs()});
  edges.resize({nx(), ny(), nnbrs()});

  build_connectivity(); 

  // The reverse inbr number.
  // For example, if the current cell is K, its 0th neighbor is L.
  // Then, for cell L, K is its 2th neighbor.
  rinbr_(0) = 2; 
  rinbr_(1) = 3; 
  rinbr_(2) = 0; 
  rinbr_(3) = 1; 
}

void Mesh::build_connectivity() {

  int inbr; 

  Point A, B; 
  Vector2 dr; 

  for (std::size_t i=0; i<nx(); ++i) {
    for (std::size_t j=0; j<ny(); ++j) {
      // nbr 0
      inbr = 0; 
      nbr_inds(i,j,inbr).i = i-1;
      nbr_inds(i,j,inbr).j = j; 

      A = {xO() + i*dx(), yO() + (j+1)*dy()}; 
      B = {xO() + i*dx(), yO() + j*dy()}; 

      edges(i,j,inbr).set_vs_dir({A,B}, XNEG); 

      // nbr 1
      inbr = 1;
      nbr_inds(i,j,inbr).i = i;
      nbr_inds(i,j,inbr).j = j+1;

      B = A;
      A = {xO() + (i+1)*dx(), yO() + (j+1)*dy()}; 

      edges(i,j,inbr).set_vs_dir({A,B}, YPOS); 

      // nbr 2
      inbr = 2; 
      nbr_inds(i,j,inbr).i = i+1;
      nbr_inds(i,j,inbr).j = j;

      B = A;
      A = {xO() + (i+1)*dx(), yO() + j*dy()}; 

      edges(i,j,inbr).set_vs_dir({A,B}, XPOS); 

      // nbr 3
      inbr = 3;
      nbr_inds(i,j,inbr).i = i;
      nbr_inds(i,j,inbr).j = j-1;

      B = A;
      A = {xO() + i*dx(), yO() + j*dy()}; 

      edges(i,j,inbr).set_vs_dir({A,B}, YNEG); 
    }
  }
}