#include "naiveMatrix.h"
#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

namespace {

    class TransposeTest : public ::testing::Test {

    protected:
        int MAX_ROWS = 100;
        int MAX_COLS = 100;
        int MAX_ELEM = 10000;
        mat_size_t n_rows;
        mat_size_t n_cols;

        std::default_random_engine generator;
        std::uniform_int_distribution<> udistC;
        std::uniform_int_distribution<> udistR;
        std::uniform_int_distribution<> udistE;

        TransposeTest() {
            srand(451);
            udistR = std::uniform_int_distribution<>(1, MAX_ROWS);
            udistC = std::uniform_int_distribution<>(1, MAX_COLS);
            udistE = std::uniform_int_distribution<>(1, MAX_ELEM);
            n_rows = static_cast<mat_size_t>(udistR(generator));
            n_cols = static_cast<mat_size_t>(udistC(generator));
        }
        template <typename T>
        Matrix<T> randomMatrix() {
            Matrix<T> m = Matrix<T>(n_rows, n_cols);
            for (int i = 0; i < n_rows; ++i)
                for (int j = 0; j < n_cols; ++j)
                    m[i][j] = udistE(generator);
            return m;
        }

        template <typename T>
        Matrix<T> randomSymmetricMatrix() {
            mat_size_t dimn = std::min(n_cols, n_rows);
            Matrix<T> m = Matrix<T>(dimn, dimn);
            for (int i = 0; i < dimn; ++i) {
                m[i][i] = static_cast<T>(udistE(generator));
                for (int j = i + 1; j < dimn; ++j) {
                    T elem = static_cast<T>(udistE(generator));
                    m[i][j] = elem;
                    m[j][i] = elem;
                }
            }
            return m;
        }
    };

    TEST_F(TransposeTest, Empty_Throws_Empty_Matrix_Ex) {
        Matrix<int> m;
        EXPECT_THROW(m.transpose(), Matrix<int>::empty_matrix);
    }

    TEST_F(TransposeTest, One_x_One_Mat_Transpose_Is_Itself) {
        Matrix<int> m(1, 1);
        Matrix<int> mT = m.transpose();
        EXPECT_EQ(m, mT);
    }

    TEST_F(TransposeTest, Transpose_Swaps_Dimensions) {
        Matrix<int> m = randomMatrix<int>();
        Matrix<int> mT = m.transpose();
        EXPECT_EQ(m.shape(0), mT.shape(1));
        EXPECT_EQ(m.shape(1), mT.shape(0));
    }

    TEST_F(TransposeTest, Sym_Matrix_Transpose_Is_Itself) {
        Matrix<int> m = randomSymmetricMatrix<int>();
        EXPECT_EQ(m, m.transpose());
    }

    TEST_F(TransposeTest, Transpose_Implementation_Equals_Naive) {
        Matrix<int> mat  = randomMatrix<int>();
        NaiveMatrix<int> naive = NaiveMatrix<int>(mat);
        EXPECT_EQ(naive.transpose(), mat.transpose());
    }
}