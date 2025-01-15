/*
 * File:        Edge.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "common.h"

enum Direction {XPOS, XNEG, YPOS, YNEG}; 

typedef std::array<Point, 2> Edge_Vertices;

class Edge{
  public:
    static const int NT = 2; // number of points for each edge. 
    const Point& v(int i) const { return vs_[i]; } // vertices
    int v_size() const { return vs_.size(); }                                                   
    double length() const { return length_; }
    const Vector2& n() const { return n_; } // normal vector, depending on the edge direction
    Direction dir() const { return dir_; }


    void set_vs_dir(const Edge_Vertices& vs, Direction dir) {
      vs_ = vs; 
      length_ = (vs_[1] - vs_[0]).norm(); 
      dir_ = dir; 

      switch (dir) {
        case XPOS: // Edge perpendicular to positive x-axis
          n_ = {1, 0};
          break;

        case XNEG:
          n_ = {-1, 0};
          break; 

        case YPOS: // Edge perpendicular to y-axis
          n_ = {0, 1};
          break;

        case YNEG:
          n_ = {0, -1};
          break;

        default: // incorrect dir
          n_ = {0, 0}; 
      }
    }

  private:

    Edge_Vertices vs_; 
    double length_; 
    Vector2 n_; // the unit normal vector
    Direction dir_; 
};

#endif /* EDGE_H */

