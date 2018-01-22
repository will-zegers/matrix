#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <vector>
#include <cstdint>
#include <sstream>
#include <tuple>
#include <iostream>

#define STEP 2
#define BLOCK_SZ 32

typedef uint32_t mat_size_t;
typedef std::tuple<mat_size_t, mat_size_t> shape_t;

template <typename T>
class Matrix {
public:
    Matrix() : n_rows(0), n_cols(0), elements(std::vector<T> ()) {}
    Matrix(shape_t shape) :
            n_rows(std::get<0>(shape)),
            n_cols(std::get<1>(shape)),
            elements(std::vector<T> (n_rows * n_cols)) {}
    Matrix(shape_t shape, std::vector<T> elements) :
            n_rows(std::get<0>(shape)),
            n_cols(std::get<1>(shape)),
            elements(elements) {}

    T& operator()(mat_size_t i, mat_size_t j) {
        return elements[i * n_cols + j];
    }
    T& operator()(mat_size_t i, mat_size_t j) const {
        return elements[i * n_cols + j];
    }

    typename std::vector<T>::const_iterator cbegin() const {
        return elements.cbegin();
    }

    typename std::vector<T>::const_iterator cend() const {
        return elements.cend();
    }

    bool empty() {
        return elements.empty();
    }

    virtual Matrix<T> transpose() {
        if (this->empty())
            throw empty_matrix();

        Matrix<T> res = Matrix<T>(std::make_pair(n_cols, n_rows));
        for (mat_size_t i = 0; i < n_rows; ++i)
            for (mat_size_t j = 0; j < n_cols; ++j)
                res(j, i) = elements[n_cols*i + j];
        return res;
    }

    virtual Matrix<T> operator*(Matrix<T>& mat) {
        if (this->empty() || mat.empty())
            throw empty_matrix();
        if (this->n_cols != mat.shape(0)) {
            throw size_mismatch();
        }

        Matrix<T> res = Matrix<T>(std::make_pair(this->n_rows, mat.shape(1)));
        mat_size_t ii;
        for (ii = 0; ii + BLOCK_SZ < this->n_rows; ii += BLOCK_SZ) {
            mat_size_t j;
            for (j = 0; j < mat.shape(1) - (STEP - 1); j += STEP) {
                mat_size_t i;
                for (i = ii; i < ii + BLOCK_SZ; i += STEP) {
                    this->blockDot2x2(i, j, mat, res);
                }  // i
                for (; i < ii + BLOCK_SZ; ++i) {
                    this->blockDot1x2(i, j, mat, res);
                }  // i
            }  // j
            for (; j < mat.shape(1); ++j) {
                mat_size_t i;
                for (i = ii; i < ii + BLOCK_SZ; i += STEP) {
                    this->blockDot2x1(i, j, mat, res);
                }  // i
                for (; i < ii + BLOCK_SZ; ++i) {
                    this->dot(i, j, mat, res);
                }  //
            }  // j
        }  // ii
        mat_size_t j;
        for (j = 0; j < mat.shape(1) - (STEP - 1); j += STEP) {
            mat_size_t i;
            for (i = ii; i < this->n_rows - (STEP - 1); i += STEP) {
                this->blockDot2x2(i, j, mat, res);
            }  // i
            for (; i < this->n_rows; ++i) {
                this->blockDot1x2(i, j, mat, res);
            }  // i
        }  // j
        for (; j < mat.shape(1); ++j) {
            mat_size_t i;
            for (i = ii; i < this->n_rows - (STEP - 1); i += STEP) {
                this->blockDot2x1(i, j, mat, res);
            }  // i
            for (; i < this->n_rows; ++i) {
                this->dot(i, j, mat, res);
            }  //
        }  // j

        return res;
    }

    bool operator==(const Matrix<T>& mat) const {
        if (n_rows != mat.shape(0) || n_cols != mat.shape(1))
            return false;

        typename std::vector<T>::const_iterator it1 = this->cbegin();
        typename std::vector<T>::const_iterator it2 = mat.cbegin();
        while (it1 != this->cend()) {
            if (*it1 != *it2)
                return false;
            ++it1; ++it2;
        }
        return true;
    }

    mat_size_t shape(mat_size_t n) const{
        return (n == 0) ? n_rows :
               (n == 1) ? n_cols :
               throw bad_shape();
    }

    std::string to_string() {
        std::stringstream ss;
        for (int i = 0; i < n_rows; ++i) {
            for (int j = 0; j < n_cols; ++j)
                ss << elements[i * n_cols + j] << "\t";
            ss << "\n";
        }

        return ss.str();
    }

    struct empty_matrix : public std::exception {
        const char* what() const throw() final {
            return "Cannot perform operation on an empty matrix";
        }
    };

    struct size_mismatch : public std::exception {
        const char* what() const throw() final {
            return "Matrix dimensions do not match";
        }
    };

    struct bad_shape : public std::exception {
        const char* what() const throw() final {
            return "Requested dimension does not exist";
        }
    };

protected:
    mat_size_t n_rows, n_cols;
    std::vector<T> elements;

    inline void blockDot2x2(mat_size_t i, mat_size_t j, Matrix<T>& other, Matrix<T>& res) {
        register T acc00 = 0, acc01 = 0, acc10 = 0, acc11 = 0;
        for (mat_size_t k = 0; k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01 += elements[n_cols * (i + 0) + k] * other(k, j + 1);
            acc10 += elements[n_cols * (i + 1) + k] * other(k, j + 0);
            acc11 += elements[n_cols * (i + 1) + k] * other(k, j + 1);
        }  // k
        res(i + 0, j + 0) = acc00;
        res(i + 0, j + 1) = acc01;
        res(i + 1, j + 0) = acc10;
        res(i + 1, j + 1) = acc11;
    }

    inline void blockDot1x2(mat_size_t i, mat_size_t j, Matrix<T>& other, Matrix<T>& res) {
        register T acc00 = 0, acc01 = 0;
        for (mat_size_t k = 0; k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01 += elements[n_cols * (i + 0) + k] * other(k, j + 1);
        }  // k
        res(i + 0, j + 0) = acc00;
        res(i + 0, j + 1) = acc01;
    }

    inline void blockDot2x1(mat_size_t i, mat_size_t j, Matrix<T>& other, Matrix<T>& res) {
        register T acc00 = 0, acc10 = 0;
        for (mat_size_t k = 0; k < this->n_cols; ++k) {
            acc00 += elements[n_cols*(i + 0) + k] * other(k, j + 0);
            acc10 += elements[n_cols*(i + 1) + k] * other(k, j + 0);
        }  // k
        res(i + 0, j + 0) = acc00;
        res(i + 1, j + 0) = acc10;
    }

    inline void dot(mat_size_t i, mat_size_t j, Matrix<T>& other, Matrix<T>& res) {
        register T acc00 = 0;
        for (mat_size_t k = 0; k < this->n_cols; ++k) {
            acc00 += elements[n_cols*(i + 0) + k] * other(k, j + 0);
        }  // k
        res(i + 0, j + 0) = acc00;
    }
};



#endif //INCLUDE_MATRIX_H
