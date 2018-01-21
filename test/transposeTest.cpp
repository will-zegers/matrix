#include "naiveMatrix.h"
#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

namespace {

    class TransposeTest : public ::testing::Test {

    protected:
        typedef int data_t;

        const int MAX_ROWS = 100;
        const int MAX_COLS = 100;
        const int MAX_DATA = 10000;
        mat_size_t dim1, dim2;

        std::default_random_engine generator;
        std::uniform_int_distribution<> uniformDim1;
        std::uniform_int_distribution<> uniformDim2;
        std::uniform_int_distribution<> uniformData;

        TransposeTest() {
            uniformDim1 = std::uniform_int_distribution<>(1, MAX_ROWS);
            uniformDim2 = std::uniform_int_distribution<>(1, MAX_COLS);
            uniformData = std::uniform_int_distribution<>(1, MAX_DATA);
            dim1 = static_cast<mat_size_t>(uniformDim1(generator));
            dim2 = static_cast<mat_size_t>(uniformDim2(generator));
        }

        Matrix<data_t> randomMatrix() {
            Matrix<data_t> m = Matrix<data_t>(std::make_pair(dim1, dim2));
            for (int i = 0; i < dim1; ++i)
                for (int j = 0; j < dim2; ++j)
                    m(i, j) = uniformData(generator);
            return m;
        }

        Matrix<data_t> randomSymmetricMatrix() {
            mat_size_t dimn = std::min(dim2, dim1);
            Matrix<data_t> m = Matrix<data_t>(std::make_pair(dimn, dimn));
            for (mat_size_t i = 0; i < dimn; ++i) {
                m(i, i) = static_cast<data_t>(uniformData(generator));
                for (mat_size_t j = i + 1; j < dimn; ++j) {
                    auto elem = static_cast<data_t>(uniformData(generator));
                    m(i, j) = elem;
                    m(j, i) = elem;
                }
            }
            return m;
        }
    };

    TEST_F(TransposeTest, Empty_Throws_Empty_Matrix_Ex) {
        Matrix<data_t> m;
        EXPECT_THROW(m.transpose(), Matrix<data_t>::empty_matrix);
    }

    TEST_F(TransposeTest, One_x_One_Mat_Transpose_Is_Itself) {
        Matrix<data_t> m = Matrix<data_t>(std::make_pair(1, 1));
        Matrix<data_t> mT = m.transpose();
        EXPECT_EQ(m, mT);
    }

    TEST_F(TransposeTest, Transpose_Swaps_Dimensions) {
        Matrix<data_t> m = randomMatrix();
        Matrix<data_t> mT = m.transpose();
        EXPECT_EQ(m.shape(0), mT.shape(1));
        EXPECT_EQ(m.shape(1), mT.shape(0));
    }

    TEST_F(TransposeTest, Sym_Matrix_Transpose_Is_Itself) {
        Matrix<data_t> m = randomSymmetricMatrix();
        EXPECT_EQ(m, m.transpose());
    }

    TEST_F(TransposeTest, Transpose_Implementation_Equals_Naive) {
        Matrix<data_t> mat = randomMatrix();
        NaiveMatrix<data_t> naive(mat);
        EXPECT_EQ(naive.transpose(), mat.transpose());
    }
}