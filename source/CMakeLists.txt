cmake_minimum_required(VERSION 2.8)

project(ppfv)

set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)
add_executable(app Parameters.hpp utils.hpp common.hpp Grid.h D.h D.cc FVMSolver.hpp FVMSolver.cc BoundaryConditions.h BoundaryConditions.cc main.cc)

if(OpenMP_CXX_FOUND)
    target_link_libraries(app OpenMP::OpenMP_CXX)
endif()
