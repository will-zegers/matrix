#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

namespace {

    class MatrixInstantiation : public ::testing::Test {
    protected:
        int MAX_ROWS = 100;
        int MAX_COLS = 100;
        mat_size_t n_rows;
        mat_size_t n_cols;

        std::default_random_engine generator;
        std::uniform_int_distribution<> udistC;
        std::uniform_int_distribution<> udistR;

        MatrixInstantiation() {
            srand(451);
            udistR = std::uniform_int_distribution<>(1, MAX_ROWS);
            udistC = std::uniform_int_distribution<>(1, MAX_COLS);
            n_rows = static_cast<mat_size_t>(udistR(generator));
            n_cols = static_cast<mat_size_t>(udistC(generator));
        }
    };

    TEST_F(MatrixInstantiation, Empty) {
        Matrix<int> m;
        EXPECT_EQ(m.shape(0), 0);
        EXPECT_EQ(m.shape(1), 0);
    }

    TEST_F(MatrixInstantiation, Empty_Is_Empty) {
        Matrix<int> m;
        EXPECT_TRUE(m.empty());
    }

    TEST_F(MatrixInstantiation, Sized) {
        Matrix<int> m(n_rows, n_cols);
        ASSERT_FALSE(m.empty());
        EXPECT_EQ(m.shape(0), n_rows);
        EXPECT_EQ(m.shape(1), n_cols);
    }

    TEST_F(MatrixInstantiation, Sized_Is_All_Zeros) {
        Matrix<int> m(n_rows, n_cols);
        for (int i = 0; i < round(sqrt(n_rows*n_cols)); ++i) {
            auto j = udistR(generator) % n_rows;
            auto k = udistC(generator) % n_cols;
            EXPECT_EQ(m[j][k], 0);
        }
    }

    TEST_F(MatrixInstantiation, From_Vec) {
        std::vector<std::vector<int> > vec =
                std::vector<std::vector<int> >(n_rows, std::vector<int>(n_cols));
        Matrix<int> m(vec);
        ASSERT_FALSE(m.empty());
        EXPECT_EQ(m.shape(0), n_rows);
        EXPECT_EQ(m.shape(1), n_cols);
    }

    TEST_F(MatrixInstantiation, From_Vec_Elems_in_Row_Major) {
        std::vector<std::vector<int> > vec =
                std::vector<std::vector<int> >(n_rows, std::vector<int>(n_cols));
        int k = 0;
        for (int i = 0; i < n_rows; ++i)
            for (int j = 0; j < n_cols; ++j)
                vec[i][j] = ++k;
        Matrix<int> m(vec);
        for (int i = 0; i < n_rows - 1; ++i) {
            EXPECT_EQ(m[i][0]+1, m[i][1]);
            EXPECT_EQ(m[i][0]+n_cols, m[i+1][0]);
            EXPECT_EQ(m[i][n_cols-2]+1, m[i][n_cols-1]);
            EXPECT_EQ(m[i][n_cols-1]+n_cols, m[i+1][n_cols-1]);
        }
    }
}