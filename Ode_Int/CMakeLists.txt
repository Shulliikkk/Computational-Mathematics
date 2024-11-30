cmake_minimum_required(VERSION 3.14)
project(ODEINT)

set(CMAKE_CXX_STANDART 17)
set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_BUILD_TYPE Release)

include_directories(src)

add_executable(odeint tests/test.cpp)
