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
#include "Mesh.h"
#include "D.h"
#include "BCs.h"
#include "Parameters.h"
#include "Solver.h"
#include <ctime>

int main() {

  // Create grid object
  Mesh m(nx, ny, dx, dy, dt);

  // Create diffusion coefficients object
  D diffusion(m);
  diffusion.constructD(0.0);

  // TODO BoundaryConditions modularization
  BCs boundary(0, 0.0);

  Solver solver(m, diffusion, boundary);

  string path;

  // The timer
  clock_t start, end;
  double cpu_time;
  start = clock();

  // Time loop for solving
  for (int k = 0; k < steps; ++k) {

    // Solve using FVM solver
    solver.update();

    if(k % printstep == printstep - 1){

      std::cout << "Time step: " << k << std::endl;
      path = "./output/SMPPFV/smppfv" + std::to_string(int((k + 1) / printstep));

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

