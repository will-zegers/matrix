#ifndef INCLUDE_NAIVEMATRIX_H
#define INCLUDE_NAIVEMATRIX_H

#include "matrix.h"

template <typename T>
class NaiveMatrix : public Matrix<T> {
public:

    explicit NaiveMatrix(std::vector<std::vector<int> > _elements) : Matrix<T>(_elements) {}
    explicit NaiveMatrix(Matrix<T> mat) : Matrix<T>(mat.all()) {}

    Matrix<T> transpose() const {
        if (this->empty())
            throw empty_matrix();

        Matrix<T> mT(this->_shape[1], this->_shape[0]);
        for (mat_size_t i = 0; i < this->_shape[1]; ++i)
            for (mat_size_t j = 0; j < this->_shape[0]; ++j)
                mT[i][j] = this->elements[j][i];

        return mT;
    }

    Matrix<T> operator*(NaiveMatrix& m) const {
        if (this->empty() || m.empty())
            throw empty_matrix();
        if (this->_shape[1] != m.shape(0))
            throw size_mismatch();

        Matrix<T> res(this->_shape[0], m.shape(1));
        for (mat_size_t i = 0; i < this->_shape[0]; ++i)
            for (mat_size_t k = 0; k < m.shape(1); ++k)
                for (mat_size_t j = 0; j < this->_shape[1]; ++j)
                    res[i][k] += this->elements[i][j] * m[j][k];
        return res;
    }

    struct empty_matrix : public std::exception {
        const char* what () const throw() {
            return "Cannot perform operation on empty matrices";
        }
    };

    struct size_mismatch : public std::exception {
        const char* what () const throw() {
            return "Incompatible matrix dimensions";
        }
    };
};

#endif //INCLUDE_NAIVEMATRIX_H
