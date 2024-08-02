/*
 * File:        Mesh.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        07/19/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Mesh.h"


void Mesh::build_connectivity() {

      int inbr; 

      Point A; 
      Point B; 

      for (std::size_t i=0; i<nx(); ++i) {
        for (std::size_t j=0; j<ny(); ++j) {
          // nbr 0
            inbr = 0; 
            nbr_inds(i,j,inbr).i = i-1;
            nbr_inds(i,j,inbr).j = j; 

            B(0) = xO_ + i*dx();  
            B(1) = yO_ + j*dy(); 

            A(0) = B(0); 
            A(1) = B(1)+dy(); 

            edges(i,j,inbr).A = A;
            edges(i,j,inbr).B = B;

          // nbr 1
            inbr = 1;
            nbr_inds(i,j,inbr).i = i;
            nbr_inds(i,j,inbr).j = j+1;
  
            B = A; 
            A(0) = B(0) + dx(); 
            A(1) = B(1); 

            edges(i,j,inbr).A = A;
            edges(i,j,inbr).B = B; 

          // nbr 2
            inbr = 2; 
            nbr_inds(i,j,inbr).i = i+1;
            nbr_inds(i,j,inbr).j = j;

            B = A; 
            A(0) = B(0);  
            A(1) = B(1) - dy();  

            edges(i,j,inbr).A = A;
            edges(i,j,inbr).B = B; 

          // nbr 3
            inbr = 3;
            nbr_inds(i,j,inbr).i = i;
            nbr_inds(i,j,inbr).j = j-1;

            B = A; 
            A(0) = B(0) - dx(); 
            A(1) = B(1); 

            edges(i,j,inbr).A = A;
            edges(i,j,inbr).B = B; 
        }
      }

    }


