cmake_minimum_required(VERSION 3.9)
project(src)

set(CMAKE_CXX_STANDARD 11)

include_directories(src)
include_directories(googletest/include googletest)

add_subdirectory(src)
add_subdirectory(tests/googletest)

add_executable(main main.cpp)

add_executable(matrixTest tests/matrixTest.cpp)
target_link_libraries(matrixTest gtest gtest_main)
