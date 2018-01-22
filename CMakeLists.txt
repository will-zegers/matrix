cmake_minimum_required(VERSION 3.9)
project(include)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)
include_directories(test/include)
include_directories(googletest/include googletest)

add_subdirectory(test/googletest)

add_executable(main main.cpp include/matrix.h)

add_executable(
        testAll
        test/include/naiveMatrix.h
        test/performanceTest.cpp)
#        test/matMulTest.cpp
#        test/transposeTest.cpp
#        test/instantiationTest.cpp)
target_link_libraries(testAll gtest gtest_main)

#SET(CMAKE_CXX_FLAGS "-O0")
#SET(CMAKE_C_FLAGS "-O0")