cmake_minimum_required(VERSION 3.9)
project(matrix)

set(CMAKE_CXX_STANDARD 11)

include_directories(googletest/include googletest)

add_subdirectory(googletest)

add_executable(matrix main.cpp matrix.cpp matrix.h main.cpp matrixTest.cpp)

add_executable(matrixTest matrixTest.cpp)
target_link_libraries(matrixTest gtest gtest_main)
