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
        typedef float data_t;

        mat_size_t dim1, dim2;
        std::chrono::duration<double, std::milli> optimDuration, naiveDuration;
        Matrix<data_t> optim1;
        NaiveMatrix<data_t> naive1;

        std::default_random_engine generator;
        std::uniform_real_distribution<> uniformData;

        PerformanceTest() {
            std::cout << std::fixed;

            generator = std::default_random_engine( (unsigned int)time(0) );
            uniformData = std::uniform_real_distribution<>(0, 1);
        }

        Matrix<data_t> randomMatrix(mat_size_t dim1, mat_size_t dim2) {
            Matrix<data_t> m = Matrix<data_t>(std::make_pair(dim1, dim2));
            for (mat_size_t i = 0; i < dim1; ++i)
                for (mat_size_t j = 0; j < dim2; ++j)
                    m(i, j) = static_cast<data_t>(uniformData(generator));
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
        dim1 = 4096, dim2 = 4096;

        std::cout << "[Info      ] " << "Testing transpose performance on size " << dim1 << " x " << dim2 << std::endl;

        optim1 = randomMatrix(dim1, dim2);
        auto start = std::chrono::high_resolution_clock::now();
        transposeWrapper(optim1);
        auto end = std::chrono::high_resolution_clock::now();
        optimDuration = end - start;

        naive1 = NaiveMatrix<data_t>(optim1);
        start = std::chrono::high_resolution_clock::now();
        transposeWrapper(naive1);
        end = std::chrono::high_resolution_clock::now();
        naiveDuration = end - start;

        std::cout << "[Naive     ] " << naiveDuration.count() << std::endl;
        std::cout << "[Optimized ] " << optimDuration.count() << std::endl;
        std::cout << "[Improvemt.] " << (naiveDuration.count() / optimDuration.count() - 1) * 100 << " %" << std::endl;

        SUCCEED();
    }

    TEST_F(PerformanceTest, MatMul) {
        dim1 = 1024, dim2 = 1024;

        std::cout << "[Info      ] " << "Testing matmul performance on size " << dim1 << " x " << dim2 << std::endl;
        optim1 = randomMatrix(dim1, dim2);
        auto start = std::chrono::high_resolution_clock::now();
        matMulWrapper(optim1);
        auto end = std::chrono::high_resolution_clock::now();
        optimDuration = end - start;

        naive1 = NaiveMatrix<data_t>(optim1);
        start = std::chrono::high_resolution_clock::now();
        matMulWrapper(naive1);
        end = std::chrono::high_resolution_clock::now();
        naiveDuration = end - start;

        std::cout << "[Naive     ] " << naiveDuration.count() << std::endl;
        std::cout << "[Optimized ] " << optimDuration.count() << std::endl;
        std::cout << "[Improvemt.] " << (naiveDuration.count() / optimDuration.count() - 1) * 100 << " %" << std::endl;

        SUCCEED();
    }
}