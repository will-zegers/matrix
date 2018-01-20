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
        std::uniform_int_distribution<> uniformDim1;
        std::uniform_int_distribution<> uniformDim2;
        std::uniform_int_distribution<> uniformDim3;
        std::uniform_int_distribution<> uniformData;

        MatMulTest() {
            uniformDim1 = std::uniform_int_distribution<>(1, MAX_DIM1);
            uniformDim2 = std::uniform_int_distribution<>(1, MAX_DIM2);
            uniformDim3 = std::uniform_int_distribution<>(1, MAX_DIM3);
            uniformData = std::uniform_int_distribution<>(1, MAX_ELEM);
            dim1 = static_cast<mat_size_t>(uniformDim1(generator));
            dim2 = static_cast<mat_size_t>(uniformDim2(generator));
            dim3 = static_cast<mat_size_t>(uniformDim3(generator));
        }

        void reroll() {
            dim1 = static_cast<mat_size_t>(uniformDim2(generator));
            dim2 = static_cast<mat_size_t>(uniformDim1(generator));
            dim3 = static_cast<mat_size_t>(uniformDim3(generator));
        }

        Matrix<data_t> randomMatrix(mat_size_t n_rows, mat_size_t n_cols) {
            Matrix<data_t> m(n_rows, n_cols);
            for (mat_size_t i = 0; i < n_rows; ++i)
                for (mat_size_t j = 0; j < n_cols; ++j)
                    m[i][j] = uniformData(generator);
            return m;
        }
    };

    TEST_F(MatMulTest, Empty_Times_Empty_Throws_Empty_Mat_Ex) {
        Matrix<data_t> m1, m2;
        EXPECT_THROW(m1 * m2, Matrix<data_t>::empty_matrix);
    }

    TEST_F(MatMulTest, Non_empty_Times_Empty_Throws_Empty_Mat_Ex) {
        Matrix<data_t> m1 = randomMatrix(dim1, dim2);
        Matrix<data_t> m2;
        EXPECT_THROW(m1 * m2, Matrix<data_t>::empty_matrix);
    }

    TEST_F(MatMulTest, Mats_with_Different_Dimensions_Throws_Size_Ex) {
        Matrix<data_t> m1 = randomMatrix(dim1, dim2);
        Matrix<data_t> m2 = randomMatrix(dim3, dim1);
        while (m1.shape(1) == m2.shape(0)) {
            reroll();
            m2 = randomMatrix(dim1, dim2);
        }
        EXPECT_THROW(m1 * m2, Matrix<data_t>::size_mismatch);
    }

    TEST_F(MatMulTest, MatMul_Implementation_Equals_Naive) {
        Matrix<data_t> mat1 = randomMatrix(dim1, dim2);
        Matrix<data_t> mat2 = randomMatrix(dim2, dim3);
        NaiveMatrix<data_t> naive1(mat1);
        NaiveMatrix<data_t> naive2(mat2);
        EXPECT_EQ(mat1 * mat2, naive1 * naive2);
    }
}