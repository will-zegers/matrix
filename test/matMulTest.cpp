#include "naiveMatrix.h"
#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

namespace {

    class MatMulTest : public ::testing::Test {

    protected:
        typedef long data_t;

        const int MAX_DIM1 = 100;
        const int MAX_DIM2 = 100;
        const int MAX_DIM3 = 100;
        const int MAX_ELEM = 10000;

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
            Matrix<T> m(n_rows, n_cols);
            for (mat_size_t i = 0; i < n_rows; ++i)
                for (mat_size_t j = 0; j < n_cols; ++j)
                    m[i][j] = udistE(generator);
            return m;
        }
    };

    TEST_F(MatMulTest, Empty_Times_Empty_Throws_Empty_Mat_Ex) {
        Matrix<data_t> m1, m2;
        EXPECT_THROW(m1 * m2, Matrix<data_t>::empty_matrix);
    }

    TEST_F(MatMulTest, Non_empty_Times_Empty_Throws_Empty_Mat_Ex) {
        Matrix<data_t> m1 = randomMatrix<data_t>(dim1, dim2);
        Matrix<data_t> m2;
        EXPECT_THROW(m1 * m2, Matrix<data_t>::empty_matrix);
    }

    TEST_F(MatMulTest, Mats_with_Different_Dimensions_Throws_Size_Ex) {
        Matrix<data_t> m1 = randomMatrix<data_t>(dim1, dim2);
        Matrix<data_t> m2 = randomMatrix<data_t>(dim3, dim1);
        while (m1.shape(1) == m2.shape(0)) {
            reroll();
            m2 = randomMatrix<data_t>(dim1, dim2);
        }
        EXPECT_THROW(m1 * m2, Matrix<data_t>::size_mismatch);
    }

    TEST_F(MatMulTest, MatMul_Implementation_Equals_Naive) {
        Matrix<data_t> mat1 = randomMatrix<data_t>(dim1, dim2);
        Matrix<data_t> mat2 = randomMatrix<data_t>(dim2, dim3);
        NaiveMatrix<data_t> naive1(mat1);
        NaiveMatrix<data_t> naive2(mat2);
        EXPECT_EQ(mat1 * mat2, naive1 * naive2);
    }
}