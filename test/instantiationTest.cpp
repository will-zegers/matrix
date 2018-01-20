#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

namespace {

    class MatrixInstantiation : public ::testing::Test {
    protected:
        typedef int data_t;
        typedef std::vector<std::vector<data_t> > vector2D;

        int MAX_ROWS = 100;
        int MAX_COLS = 100;
        mat_size_t dim1, dim2;

        std::default_random_engine generator;
        std::uniform_int_distribution<> uniformDim1;
        std::uniform_int_distribution<> uniformDim2;

        MatrixInstantiation() {
            uniformDim1 = std::uniform_int_distribution<>(1, MAX_ROWS);
            uniformDim2 = std::uniform_int_distribution<>(1, MAX_COLS);
            dim1 = static_cast<mat_size_t>(uniformDim1(generator));
            dim2 = static_cast<mat_size_t>(uniformDim2(generator));
        }
    };

    TEST_F(MatrixInstantiation, Empty) {
        Matrix<data_t> m;
        EXPECT_EQ(m.shape(0), 0);
        EXPECT_EQ(m.shape(1), 0);
    }

    TEST_F(MatrixInstantiation, Empty_Is_Empty) {
        Matrix<data_t> m;
        EXPECT_TRUE(m.empty());
    }

    TEST_F(MatrixInstantiation, Sized) {
        Matrix<data_t> m(dim1, dim2);
        ASSERT_FALSE(m.empty());
        EXPECT_EQ(m.shape(0), dim1);
        EXPECT_EQ(m.shape(1), dim2);
    }

    TEST_F(MatrixInstantiation, Sized_Is_All_Zeros) {
        Matrix<int> m(dim1, dim2);
        for (mat_size_t i = 0; i < round(sqrt(dim1*dim2)); ++i) {
            auto j = uniformDim1(generator) % dim1;
            auto k = uniformDim2(generator) % dim2;
            EXPECT_EQ(m[j][k], 0);
        }
    }

    TEST_F(MatrixInstantiation, From_Vec) {
        vector2D vec = vector2D(dim1, std::vector<data_t>(dim2));
        Matrix<data_t> m(vec);
        ASSERT_FALSE(m.empty());
        EXPECT_EQ(m.shape(0), dim1);
        EXPECT_EQ(m.shape(1), dim2);
    }

    TEST_F(MatrixInstantiation, From_Vec_Elems_in_Row_Major) {
        vector2D vec = vector2D(dim1, std::vector<data_t>(dim2));
        data_t k = 0;
        for (mat_size_t i = 0; i < dim1; ++i)
            for (mat_size_t j = 0; j < dim2; ++j)
                vec[i][j] = ++k;
        Matrix<data_t> m(vec);
        for (mat_size_t i = 0; i < dim1 - 1; ++i) {
            EXPECT_EQ(m[i][0]+1, m[i][1]);
            EXPECT_EQ(m[i][0]+dim2, m[i+1][0]);
            EXPECT_EQ(m[i][dim2-2]+1, m[i][dim2-1]);
            EXPECT_EQ(m[i][dim2-1]+dim2, m[i+1][dim2-1]);
        }
    }
}