#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <vector>
#include <cstdint>
#include <sstream>
#include <tuple>
#include <iostream>

#define XPOSE_STEP 2
#define MATMUL_STEP 4
#define I_BLOCK_SZ 64
#define K_BLOCK_SZ 64

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

    inline T& operator()(mat_size_t i, mat_size_t j) {
        return elements[i * n_cols + j];
    }
    inline T& operator()(mat_size_t i, mat_size_t j) const {
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
        mat_size_t i;
        for (i = 0; i + XPOSE_STEP < n_rows; i += XPOSE_STEP) {
            mat_size_t j;
            for (j = 0; j + XPOSE_STEP < n_cols; j += XPOSE_STEP) {
                res(j + 0, i + 0) = elements[n_cols * (i + 0) + (j + 0)];
                res(j + 1, i + 0) = elements[n_cols * (i + 0) + (j + 1)];
                res(j + 0, i + 1) = elements[n_cols * (i + 1) + (j + 0)];
                res(j + 1, i + 1) = elements[n_cols * (i + 1) + (j + 1)];
            }  // j
            for (; j < n_cols; ++j) {
                res(j + 0, i + 0) = elements[n_cols * (i + 0) + (j + 0)];
                res(j + 0, i + 1) = elements[n_cols * (i + 1) + (j + 0)];
            }  // j
        }  // i
        for (; i < n_rows; ++i) {
            mat_size_t j;
            for (j = 0; j + XPOSE_STEP < n_cols; j += XPOSE_STEP) {
                res(j + 0, i + 0) = elements[n_cols * (i + 0) + (j + 0)];
                res(j + 1, i + 0) = elements[n_cols * (i + 0) + (j + 1)];
            }  // j
            for (; j < n_cols; ++j) {
                res(j + 0, i + 0) = elements[n_cols * (i + 0) + (j + 0)];
            }  // j
        }  // i
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
        for (ii = 0; ii + I_BLOCK_SZ < this->n_rows; ii += I_BLOCK_SZ) {
            mat_size_t kk;
            for (kk = 0; kk + K_BLOCK_SZ < this->n_cols; kk += K_BLOCK_SZ) {
                mat_size_t j;
                for (j = 0; j + MATMUL_STEP < mat.shape(1); j += MATMUL_STEP) {
                    mat_size_t i;
                    for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                        this->blockDotNxN(i, j, kk, mat, res);
                    }  // i
                    for (; i < ii + I_BLOCK_SZ; ++i) {
                        this->blockDot1xN(i, j, kk, mat, res);
                    }  // i
                }  // j
                for (; j < mat.shape(1); ++j) {
                    mat_size_t i;
                    for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                        this->blockDotNx1(i, j, kk, mat, res);
                    }  // i
                    for (; i < ii + I_BLOCK_SZ; ++i) {
                        this->dot(i, j, kk, mat, res);
                    }  // i
                }  // j
            }  // kk
            mat_size_t j;
            for (j = 0; j + MATMUL_STEP < mat.shape(1); j += MATMUL_STEP) {
                mat_size_t i;
                for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                    this->blockDotNxN(i, j, kk, mat, res);
                }  // i
                for (; i < ii + I_BLOCK_SZ; ++i) {
                    this->blockDot1xN(i, j, kk, mat, res);
                }  // i
            }  // j
            for (; j < mat.shape(1); ++j) {
                mat_size_t i;
                for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                    this->blockDotNx1(i, j, kk, mat, res);
                }  // i
                for (; i < ii + I_BLOCK_SZ; ++i) {
                    this->dot(i, j, kk, mat, res);
                }  // i
            }  // j
        }  // ii
        mat_size_t kk;
        for (kk = 0; kk + K_BLOCK_SZ < this->n_cols; kk += K_BLOCK_SZ) {
            mat_size_t j;
            for (j = 0; j + MATMUL_STEP < mat.shape(1); j += MATMUL_STEP) {
                mat_size_t i;
                for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                    this->blockDotNxN(i, j, kk, mat, res);
                }  // i
                for (; i < this->n_rows; ++i) {
                    this->blockDot1xN(i, j, kk, mat, res);
                }  // i
            }  // j
            for (; j < mat.shape(1); ++j) {
                mat_size_t i;
                for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                    this->blockDotNx1(i, j, kk, mat, res);
                }  // i
                for (; i < this->n_rows; ++i) {
                    this->dot(i, j, kk, mat, res);
                }  // i
            }  // j
        }  // kk
        mat_size_t j;
        for (j = 0; j + MATMUL_STEP < mat.shape(1); j += MATMUL_STEP) {
            mat_size_t i;
            for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                this->blockDotNxN(i, j, kk, mat, res);
            }  // i
            for (; i < this->n_rows; ++i) {
                this->blockDot1xN(i, j, kk, mat, res);
            }  // i
        }  // j
        for (; j < mat.shape(1); ++j) {
            mat_size_t i;
            for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                this->blockDotNx1(i, j, kk,  mat, res);
            }  // i
            for (; i < this->n_rows; ++i) {
                this->dot(i, j, kk, mat, res);
            }  // i
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

    inline void blockDotNxN(mat_size_t i, mat_size_t j, mat_size_t kk, Matrix<T> &other, Matrix<T> &res) {
        register T acc00 = 0, acc10 = 0, acc20 = 0, acc30 = 0, acc01 = 0, acc11 = 0, acc21 = 0, acc31 = 0;
        register T acc02 = 0, acc12 = 0, acc22 = 0, acc32 = 0, acc03 = 0, acc13 = 0, acc23 = 0, acc33 = 0;

        if (kk == 0) {
           acc00 = 0; acc10 = 0; acc20 = 0; acc30 = 0; acc01 = 0; acc11 = 0; acc21 = 0; acc31 = 0;
           acc02 = 0; acc12 = 0; acc22 = 0; acc32 = 0; acc03 = 0; acc13 = 0; acc23 = 0; acc33 = 0;
        } else {
            acc00 = res(i + 0, j + 0); acc20 = res(i + 2, j + 0);
            acc01 = res(i + 0, j + 1); acc21 = res(i + 2, j + 1);
            acc02 = res(i + 0, j + 2); acc22 = res(i + 2, j + 2);
            acc03 = res(i + 0, j + 3); acc23 = res(i + 2, j + 3);
            acc10 = res(i + 1, j + 0); acc30 = res(i + 3, j + 0);
            acc11 = res(i + 1, j + 1); acc31 = res(i + 3, j + 1);
            acc12 = res(i + 1, j + 2); acc32 = res(i + 3, j + 2);
            acc13 = res(i + 1, j + 3); acc33 = res(i + 3, j + 3);
        }

        for (mat_size_t k = kk; k < kk + K_BLOCK_SZ && k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01 += elements[n_cols * (i + 0) + k] * other(k, j + 1);
            acc02 += elements[n_cols * (i + 0) + k] * other(k, j + 2);
            acc03 += elements[n_cols * (i + 0) + k] * other(k, j + 3);

            acc10 += elements[n_cols * (i + 1) + k] * other(k, j + 0);
            acc11 += elements[n_cols * (i + 1) + k] * other(k, j + 1);
            acc12 += elements[n_cols * (i + 1) + k] * other(k, j + 2);
            acc13 += elements[n_cols * (i + 1) + k] * other(k, j + 3);

            acc20 += elements[n_cols * (i + 2) + k] * other(k, j + 0);
            acc21 += elements[n_cols * (i + 2) + k] * other(k, j + 1);
            acc22 += elements[n_cols * (i + 2) + k] * other(k, j + 2);
            acc23 += elements[n_cols * (i + 2) + k] * other(k, j + 3);

            acc30 += elements[n_cols * (i + 3) + k] * other(k, j + 0);
            acc31 += elements[n_cols * (i + 3) + k] * other(k, j + 1);
            acc32 += elements[n_cols * (i + 3) + k] * other(k, j + 2);
            acc33 += elements[n_cols * (i + 3) + k] * other(k, j + 3);
        }  // k
        res(i + 0, j + 0) = acc00;
        res(i + 0, j + 1) = acc01;
        res(i + 0, j + 2) = acc02;
        res(i + 0, j + 3) = acc03;

        res(i + 1, j + 0) = acc10;
        res(i + 1, j + 1) = acc11;
        res(i + 1, j + 2) = acc12;
        res(i + 1, j + 3) = acc13;

        res(i + 2, j + 0) = acc20;
        res(i + 2, j + 1) = acc21;
        res(i + 2, j + 2) = acc22;
        res(i + 2, j + 3) = acc23;

        res(i + 3, j + 0) = acc30;
        res(i + 3, j + 1) = acc31;
        res(i + 3, j + 2) = acc32;
        res(i + 3, j + 3) = acc33;
    }

    inline void blockDotNx1(mat_size_t i, mat_size_t j, mat_size_t kk, Matrix<T> &other, Matrix<T> &res) {
        register T acc00 = 0, acc10 = 0, acc20 = 0, acc30 = 0;

        if (kk == 0) {
            acc00 = 0; acc10 = 0; acc20 = 0; acc30 = 0;
        } else {
            acc00 = res(i + 0, j + 0);
            acc10 = res(i + 1, j + 0);
            acc20 = res(i + 2, j + 0);
            acc30 = res(i + 3, j + 0);
        }
        for (mat_size_t k = kk; k < kk + K_BLOCK_SZ && k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc10 += elements[n_cols * (i + 1) + k] * other(k, j + 0);
            acc20 += elements[n_cols * (i + 2) + k] * other(k, j + 0);
            acc30 += elements[n_cols * (i + 3) + k] * other(k, j + 0);
        }

        res(i + 0, j + 0) = acc00;
        res(i + 1, j + 0) = acc10;
        res(i + 2, j + 0) = acc20;
        res(i + 3, j + 0) = acc30;
    }

    inline void blockDot1xN(mat_size_t i, mat_size_t j, mat_size_t kk, Matrix<T> &other, Matrix<T> &res) {
        register T acc00 = 0, acc01 = 0, acc02 = 0, acc03 = 0;

        if (kk == 0) {
            acc00 = 0;acc01 = 0; acc02 = 0;acc03 = 0;
        } else {
            acc00 = res(i + 0, j + 0);
            acc01 = res(i + 0, j + 1);
            acc02 = res(i + 0, j + 2);
            acc03 = res(i + 0, j + 3);
        }
        for (mat_size_t k = kk; k < kk + K_BLOCK_SZ && k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01 += elements[n_cols * (i + 0) + k] * other(k, j + 1);
            acc02 += elements[n_cols * (i + 0) + k] * other(k, j + 2);
            acc03 += elements[n_cols * (i + 0) + k] * other(k, j + 3);
        }  // k
        res(i + 0, j + 0) = acc00;
        res(i + 0, j + 1) = acc01;
        res(i + 0, j + 2) = acc02;
        res(i + 0, j + 3) = acc03;
    }

    inline void dot(mat_size_t i, mat_size_t j, mat_size_t kk, Matrix<T>& other, Matrix<T>& res) {
        register T acc00 = 0;

        acc00 = (kk == 0) ? 0 : res(i + 0, j + 0);

        for (mat_size_t k = kk; k < kk + K_BLOCK_SZ && k < this->n_cols; ++k) {
            acc00 += elements[n_cols*(i + 0) + k] * other(k, j + 0);
        }  // k
        res(i + 0, j + 0) = acc00;
    }
};



#endif //INCLUDE_MATRIX_H
