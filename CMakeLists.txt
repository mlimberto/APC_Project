cmake_minimum_required(VERSION 2.8)

project("Test_Project")

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -g")

include_directories(./include)

add_executable(out src/main.cpp
	src/test_lib.cpp
	src/lib_brando.cpp
	src/testlib_lara.cpp
)