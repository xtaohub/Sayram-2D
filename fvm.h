/*
 * File:        fvm.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        01/28/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef FVM_H_
#define FVM_H_

#include "common.h"

class Grid; 
void fvm_update(const Grid& g, const Matrix& fn, Matrix* fnp1p);

#endif /* FVM_H */

