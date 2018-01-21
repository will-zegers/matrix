#include "matrix.h"
#include <gtest/gtest.h>
#include <random>

namespace {

    class MatrixInstantiation : public ::testing::Test {
    protected:
        typedef int data_t;

        int MAX_ROWS = 100;
        int MAX_COLS = 100;
        int MAX_DATA = 10000;
        mat_size_t dim1, dim2;

        std::default_random_engine generator;
        std::uniform_int_distribution<> uniformDim1;
        std::uniform_int_distribution<> uniformDim2;
        std::uniform_int_distribution<> uniformData;

        MatrixInstantiation() {
            uniformDim1 = std::uniform_int_distribution<>(1, MAX_ROWS);
            uniformDim2 = std::uniform_int_distribution<>(1, MAX_COLS);
            uniformData = std::uniform_int_distribution<>(1, MAX_DATA);
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
        Matrix<data_t> m = Matrix<data_t>(std::make_pair(dim1, dim2));
        ASSERT_FALSE(m.empty());
        EXPECT_EQ(m.shape(0), dim1);
        EXPECT_EQ(m.shape(1), dim2);
    }

    TEST_F(MatrixInstantiation, Sized_Is_All_Zeros) {
        Matrix<int> m = Matrix<data_t>(std::make_pair(dim1, dim2));
        for (mat_size_t i = 0; i < round(sqrt(dim1*dim2)); ++i) {
            auto j = uniformDim1(generator) % dim1;
            auto k = uniformDim2(generator) % dim2;
            EXPECT_EQ(m(j, k), 0);
        }
    }

    TEST_F(MatrixInstantiation, From_Vec) {
        std::vector<data_t> vec = std::vector<data_t>(dim1 * dim2);
        Matrix<data_t> m = Matrix<data_t>(std::make_pair(dim1, dim2));
        ASSERT_FALSE(m.empty());
        EXPECT_EQ(m.shape(0), dim1);
        EXPECT_EQ(m.shape(1), dim2);
    }

    TEST_F(MatrixInstantiation, From_Vec_Elems_in_Row_Major) {
        std::vector<data_t> vec = std::vector<data_t>(dim1 * dim2);
        for (mat_size_t i = 0; i < dim1; ++i)
            for (mat_size_t j = 0; j < dim2; ++j)
                vec[dim2*i+j] = dim2*i + j;

        Matrix<data_t> m = Matrix<data_t>(std::make_pair(dim1, dim2));
        for (mat_size_t i = 0; i < dim1 - 1; ++i) {
            EXPECT_EQ(m(i, 0)+1, m(i, 1));
            EXPECT_EQ(m(i, 0)+dim2, m(i+1, 0));
            EXPECT_EQ(m(i, dim2-2)+1, m(i, dim2-1));
            EXPECT_EQ(m(i, dim2-1)+dim2, m(i+1, dim2-1));
        }
    }

    TEST_F(MatrixInstantiation, Can_Add_And_Read) {
        Matrix<data_t> m = Matrix<data_t>(std::make_pair(dim1, dim2));
        auto elem = static_cast<data_t>(uniformData(generator));
        mat_size_t i = uniformDim1(generator) % dim1;
        mat_size_t j = uniformDim2(generator) % dim2;
        m(i, j) = elem;
        EXPECT_EQ(elem, m(i, j));
    }

    TEST_F(MatrixInstantiation, Accessing_Beyond_Dims_Throw_OOB) {
        Matrix<data_t> m = Matrix<data_t>(std::make_pair(dim1, dim2));
        EXPECT_THROW(m(dim1, dim2-1), Matrix<data_t>::out_of_bounds);
        EXPECT_THROW(m(dim1-1, dim2), Matrix<data_t>::out_of_bounds);
        EXPECT_THROW(m(dim1, dim2), Matrix<data_t>::out_of_bounds);
    }
}