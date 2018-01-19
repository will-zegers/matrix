#include "matrix.h"
#include "matrix.cpp"
#include <gtest/gtest.h>

TEST(Insantiate_Matrix, Empty) {
    Matrix<int> m;
    EXPECT_EQ(m.shape(0), 0);
    EXPECT_EQ(m.shape(1), 0);
    ASSERT_EXIT((m.at(0, 0)),::testing::KilledBySignal(SIGSEGV),".*");
}

TEST(Instantiate_Matrix, Zero) {
    srand(451);
    for (size_t i = 0; i < 10; ++i) {
        auto rows = static_cast<size_t>(rand() % 1000);
        auto cols = static_cast<size_t>(rand() % 1000);
        Matrix<int> m(rows, cols);
        EXPECT_EQ(m.shape(0), rows);
        EXPECT_EQ(m.shape(1), cols);
    }
}