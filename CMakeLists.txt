cmake_minimum_required(VERSION 2.8)

project("APC_Project")

set (APC_Project_VERSION_MAJOR 1)
set (APC_Project_VERSION_MAJOR 0)

IF (APPLE)
	message("-- Setting configurations for Mac OS X")
	set(CMAKE_CXX_COMPILER "clang++")
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -g -larmadillo")
ELSE()
	message("-- Setting configurations for Linux machine")
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -g -llapack -lblas -Iarmadillo-4.550.3/include/ ")
ENDIF()

include_directories(./include)
include_directories(../APC_Build/armadillo-4.550.3/include)

add_executable(out src/my_utils.cpp 
				   src/pg_u.cpp 
				   src/pg_h.cpp 
				   src/amf.cpp
				   src/main.cpp 
				   )
