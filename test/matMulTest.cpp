#include "naiveMatrix.h"
#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

namespace {

    class MatMulTest : public ::testing::Test {

    protected:
        int MAX_DIM1 = 100;
        int MAX_DIM2 = 100;
        int MAX_DIM3 = 100;
        int MAX_ELEM = 10000;
        mat_size_t dim1, dim2, dim3;

        std::default_random_engine generator;
        std::uniform_int_distribution<> udistDim1;
        std::uniform_int_distribution<> udistDim2;
        std::uniform_int_distribution<> udistDim3;
        std::uniform_int_distribution<> udistE;

        MatMulTest() {
            srand(451);
            udistDim1 = std::uniform_int_distribution<>(1, MAX_DIM1);
            udistDim2 = std::uniform_int_distribution<>(1, MAX_DIM2);
            udistDim3 = std::uniform_int_distribution<>(1, MAX_DIM3);
            udistE = std::uniform_int_distribution<>(1, MAX_ELEM);
            dim1 = static_cast<mat_size_t>(udistDim1(generator));
            dim2 = static_cast<mat_size_t>(udistDim2(generator));
            dim3 = static_cast<mat_size_t>(udistDim3(generator));
        }

        void reroll() {
            dim1 = static_cast<mat_size_t>(udistDim2(generator));
            dim2 = static_cast<mat_size_t>(udistDim1(generator));
            dim3 = static_cast<mat_size_t>(udistDim3(generator));
        }

        template <typename T>
        Matrix<T> randomMatrix(mat_size_t n_rows, mat_size_t n_cols) {
            Matrix<T> m = Matrix<T>(n_rows, n_cols);
            for (int i = 0; i < n_rows; ++i)
                for (int j = 0; j < n_cols; ++j)
                    m[i][j] = udistE(generator);
            return m;
        }
    };

    TEST_F(MatMulTest, Empty_Times_Empty_Throws_Empty_Mat_Ex) {
        Matrix<int> m1, m2;
        EXPECT_THROW(m1 * m2, Matrix<int>::empty_matrix);
    }

    TEST_F(MatMulTest, Non_empty_Times_Empty_Throws_Empty_Mat_Ex) {
        Matrix<int> m1 = randomMatrix<int>(dim1, dim2);
        Matrix<int> m2;
        EXPECT_THROW(m1 * m2, Matrix<int>::empty_matrix);
    }

    TEST_F(MatMulTest, Mats_with_Different_Dimensions_Throws_Size_Ex) {
        Matrix<int> m1 = randomMatrix<int>(dim1, dim2);
        Matrix<int> m2 = randomMatrix<int>(dim3, dim1);
        while (m1.shape(1) == m2.shape(0)) {
            reroll();
            m2 = randomMatrix<int>(dim1, dim2);
        }
        EXPECT_THROW(m1 * m2, Matrix<int>::size_mismatch);
    }

    TEST_F(MatMulTest, MatMul_Implementation_Equals_Naive) {
        Matrix<long> mat1 = randomMatrix<long>(dim1, dim2);
        Matrix<long> mat2 = randomMatrix<long>(dim2, dim3);
        NaiveMatrix<long> naive1 = NaiveMatrix<long>(mat1);
        NaiveMatrix<long> naive2 = NaiveMatrix<long>(mat2);
        EXPECT_EQ(mat1 * mat2, naive1 * naive2);
    }
}