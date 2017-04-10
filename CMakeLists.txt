cmake_minimum_required(VERSION 2.8)
project(untitled)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp compute_shape.cpp compute_dfxi.cpp compute_dj.cpp compute_dji.cpp compute_dfx.cpp main.h compute_dfx.h set_mesh.cpp set_bc.cpp set_connectivity.cpp)
include_directories("/usr/local/Cellar/eigen/3.2.6/include/eigen3/")
add_executable(untitled ${SOURCE_FILES})

