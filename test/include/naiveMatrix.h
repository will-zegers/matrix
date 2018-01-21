#ifndef INCLUDE_NAIVEMATRIX_H
#define INCLUDE_NAIVEMATRIX_H

#include "matrix.h"

template <typename T>
class NaiveMatrix : public Matrix<T> {
public:

    NaiveMatrix() : Matrix<T>() {}
    explicit NaiveMatrix(shape_t shape, std::vector<T> elements) :
            Matrix<T>(shape, elements) {}
    explicit NaiveMatrix(Matrix<T>& mat) :
            Matrix<T>(std::make_pair(mat.shape(0), mat.shape(1)), std::vector<T>(mat.cbegin(), mat.cend())) {}

    virtual Matrix<T> transpose() {
        if (this->empty())
            throw typename Matrix<T>::empty_matrix();

        Matrix<T> res = Matrix<T>(std::make_pair(this->n_cols, this->n_rows));
        for (mat_size_t i = 0; i < this->n_rows; ++i)
            for (mat_size_t j = 0; j < this->n_cols; ++j)
                res(j, i) = this->elements[this->n_cols*i + j];
        return res;
    }

    virtual Matrix<T> operator*(Matrix<T>& mat) {
        if (this->empty() || mat.empty())
            throw typename Matrix<T>::empty_matrix();
        if (this->n_cols != mat.shape(0)) {
            throw typename Matrix<T>::size_mismatch();
        }

        Matrix<T> res(std::make_pair(this->n_rows, mat.shape(1)));
        for (mat_size_t i = 0; i < this->n_rows; ++i)
            for (mat_size_t k = 0; k < mat.shape(1); ++k)
                for (mat_size_t j = 0; j < this->n_cols; ++j)
                    res(i, k) += this->elements[this->n_cols*i + j] * mat(j, k);
        return res;
    }
};

#endif //INCLUDE_NAIVEMATRIX_H
