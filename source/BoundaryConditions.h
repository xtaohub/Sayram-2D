/*
 * File:        BoundaryConditions.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "Eigen/Dense"

class BoundaryConditions {
public:
    BoundaryConditions();

    // Apply boundary conditions to the given solution matrix f
    void applyBoundaryConditions(Eigen::MatrixXd& f);

    // Define your boundary condition functions here

private:
    // Define private members and helper functions here
};

#endif /* BOUNDARY_CONDITIONS_H */
