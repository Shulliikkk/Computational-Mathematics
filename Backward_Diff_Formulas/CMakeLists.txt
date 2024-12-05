cmake_minimum_required(VERSION 3.14)
project(BDF)

set(CMAKE_CXX_STANDART 20)
set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_BUILD_TYPE Release)

add_executable(bdf tests/test.cpp)

include_directories("src")
include_directories("third_party/eigen/")
