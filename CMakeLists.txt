cmake_minimum_required(VERSION 2.8)

project("APC_Project")

set (APC_Project_VERSION_MAJOR 1)
set (APC_Project_VERSION_MAJOR 0)

IF (APPLE)
	message("-- Setting configurations for Mac OS X")
	set(CMAKE_CXX_COMPILER "clang++")
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -g -larmadillo")
	add_definitions(-DAMFTIME )
ELSE()
	message("-- Setting configurations for Linux machine")
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -g -llapack -lblas -Iarmadillo-4.550.3/include/ ")
	include_directories(../APC_Build/armadillo-4.550.3/include)
ENDIF()

include_directories(./include)

add_executable(amf_out src/my_utils.cpp
				   src/compute_v.cpp
				   src/proj_gradient.cpp 
				   src/amf.cpp
				   main/main_tuning.cpp
				   )
