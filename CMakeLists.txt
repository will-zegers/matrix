cmake_minimum_required(VERSION 3.9)
project(include)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)
include_directories(googletest/include googletest)

add_subdirectory(test/googletest)

add_executable(main main.cpp include/matrix.h)

add_executable(matrixTest test/matrixTest.cpp)
target_link_libraries(matrixTest gtest gtest_main)
