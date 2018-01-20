#include <iostream>
#include <random>
#include <iomanip>

#include "matrix.h"
#include "naiveMatrix.h"
#include <gtest/gtest.h>
#include <chrono>

namespace {
    class PerformanceTest : public ::testing::Test {

    protected:
        typedef int data_t;

        const int MAX_ELEM = 10000;

        mat_size_t dimn;
        std::chrono::duration<double, std::milli> duration;
        Matrix<data_t> optim1;
        NaiveMatrix<data_t> naive1;

        std::default_random_engine generator;
        std::uniform_int_distribution<> uniformData;

        PerformanceTest() {
            std::cout << std::fixed;
            uniformData = std::uniform_int_distribution<>(1, MAX_ELEM);
        }

        Matrix<data_t> randomMatrix(mat_size_t dimn) {
            Matrix<data_t> m(dimn, dimn);
            for (mat_size_t i = 0; i < dimn; ++i)
                for (mat_size_t j = 0; j < dimn; ++j)
                    m[i][j] = uniformData(generator);
            return m;
        }

        Matrix<data_t> transposeWrapper(Matrix<data_t>& mat) {
            return mat.transpose();
        }

        Matrix<data_t> transposeWrapper(NaiveMatrix<data_t>& mat) {
            return mat.transpose();
        }

        Matrix<data_t> matMulWrapper(Matrix<data_t>& mat) {
            return mat * mat;
        }

        Matrix<data_t> matMulWrapper(NaiveMatrix<data_t>& mat) {
            return mat * mat;
        }
    };

    TEST_F(PerformanceTest, Transpose) {
        dimn = 4096;

        optim1 = randomMatrix(dimn);
        auto start = std::chrono::high_resolution_clock::now();
        transposeWrapper(optim1);
        auto end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "[Optimized ] " << duration.count() << std::endl;

        naive1 = NaiveMatrix<data_t>(optim1);
        start = std::chrono::high_resolution_clock::now();
        transposeWrapper(naive1);
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "[Naive     ] " << duration.count() << std::endl;

        SUCCEED();
    }

    TEST_F(PerformanceTest, MatMul) {
        dimn = 256;

        optim1 = randomMatrix(dimn);
        auto start = std::chrono::high_resolution_clock::now();
        matMulWrapper(optim1);
        auto end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "[Optimized ] " << duration.count() << std::endl;

        naive1 = NaiveMatrix<data_t>(optim1);
        start = std::chrono::high_resolution_clock::now();
        matMulWrapper(naive1);
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "[Naive     ] " << duration.count() << std::endl;

        SUCCEED();
    }
}