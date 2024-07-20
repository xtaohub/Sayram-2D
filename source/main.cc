/*
 * fin:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

// #define EIGEN_USE_THREADS

// TODO MKL
// #define EIGEN_USE_MKL_ALL
// #define EIGEN_VECTORIZE_SSE4_2

#include <iostream>
#include <cassert>
#include "Parameters.h"
#include "Mesh.h"
#include "D.h"
#include "BCs.h"
#include "Solver.h"
#include <ctime>

int main(int argc, char** argv) {

  Parameters paras(argc,argv); 

  // Create grid object
  Mesh m(paras);

  // Create diffusion coefficients object
  D diffusion(paras, m);

  // TODO BoundaryConditions modularization
  BCs boundary(paras);

  Solver solver(m, diffusion, boundary);

  string path;

  // The timer
  clock_t start, end;
  double cpu_time;
  start = clock();

  // Time loop for solving
  for (int k = 1; k <= paras.nsteps(); ++k) {

    // Solve using FVM solver
    solver.update();

    if(k % paras.save_every_step() == 0){

      path = paras.output_path() + "/" + paras.run_id() + std::to_string(int((k) / paras.save_every_step()));

      ofstream outFile(path);
      assert(outFile);

      for (double value : solver.f().reshaped()) {
        outFile << value << std::endl;
      }
      outFile.close();
    }
  }
  end = clock();
  cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used " << cpu_time << " seconds" << std::endl;

  return 0;
}

