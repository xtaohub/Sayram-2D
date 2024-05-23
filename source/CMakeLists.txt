cmake_minimum_required(VERSION 2.8)

project(ppfv)

add_executable(app Parameters.hpp utils.hpp common.hpp Grid.h D.h D.cc main.cc)