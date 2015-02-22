cmake_minimum_required(VERSION 2.8)

project("APC_Project")

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -Wall -g -larmadillo")
add_definitions(-DAMFTIME -DAMFOPENMP )

include_directories(./include)

add_executable(amf_out src/my_utils.cpp
				   src/compute_v.cpp
				   src/proj_gradient.cpp 
				   src/amf.cpp
				   main/main_solve.cpp
				   )
