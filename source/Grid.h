/*
 * File:        Grid.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef GRID_H_
#define GRID_H_

class Grid {
public:
    Grid(int nx, int ny, double dx, double dy, double dt):
      nx_(nx), ny_(ny), dx_(dx), dy_(dy), dt_(dt) {}

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double dt() const { return dt_; }

private:
    int nx_;
    int ny_;
    double dx_;
    double dy_;
    double dt_; 
};

#endif /* GRID_H */

