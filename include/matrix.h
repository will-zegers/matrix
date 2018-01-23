#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <vector>
#include <cstdint>
#include <sstream>
#include <tuple>
#include <iostream>

#define XPOSE_STEP 2
#define MATMUL_STEP 8
#define BLOCK_SZ 64

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
        for (ii = 0; ii + BLOCK_SZ < this->n_rows; ii += BLOCK_SZ) {
            mat_size_t j;
            for (j = 0; j + MATMUL_STEP < mat.shape(1); j += MATMUL_STEP) {
                mat_size_t i;
                for (i = ii; i < ii + BLOCK_SZ; i += MATMUL_STEP) {
                    this->blockDotNxN(i, j, mat, res);
                }  // i
                for (; i < ii + BLOCK_SZ; ++i) {
                    this->blockDot1xN(i, j, mat, res);
                }  // i
            }  // j
            for (; j < mat.shape(1); ++j) {
                mat_size_t i;
                for (i = ii; i < ii + BLOCK_SZ; i += MATMUL_STEP) {
                    this->blockDotNx1(i, j, mat, res);
                }  // i
                for (; i < ii + BLOCK_SZ; ++i) {
                    this->dot(i, j, mat, res);
                }  // i
            }  // j
        }  // ii
        mat_size_t j;
        for (j = 0; j + MATMUL_STEP < mat.shape(1); j += MATMUL_STEP) {
            mat_size_t i;
            for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                this->blockDotNxN(i, j, mat, res);
            }  // i
            for (; i < this->n_rows; ++i) {
                this->blockDot1xN(i, j, mat, res);
            }  // i
        }  // j
        for (; j < mat.shape(1); ++j) {
            mat_size_t i;
            for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                this->blockDotNx1(i, j, mat, res);
            }  // i
            for (; i < this->n_rows; ++i) {
                this->dot(i, j, mat, res);
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

    inline void blockDotNxN(mat_size_t i, mat_size_t j, Matrix<T> &other, Matrix<T> &res) {
        register T acc00 = 0, acc10 = 0, acc20 = 0, acc30 = 0, acc40 = 0, acc50 = 0, acc60 = 0, acc70 = 0;
        register T acc01 = 0, acc11 = 0, acc21 = 0, acc31 = 0, acc41 = 0, acc51 = 0, acc61 = 0, acc71 = 0;
        register T acc02 = 0, acc12 = 0, acc22 = 0, acc32 = 0, acc42 = 0, acc52 = 0, acc62 = 0, acc72 = 0;
        register T acc03 = 0, acc13 = 0, acc23 = 0, acc33 = 0, acc43 = 0, acc53 = 0, acc63 = 0, acc73 = 0;
        register T acc04 = 0, acc14 = 0, acc24 = 0, acc34 = 0, acc44 = 0, acc54 = 0, acc64 = 0, acc74 = 0;
        register T acc05 = 0, acc15 = 0, acc25 = 0, acc35 = 0, acc45 = 0, acc55 = 0, acc65 = 0, acc75 = 0;
        register T acc06 = 0, acc16 = 0, acc26 = 0, acc36 = 0, acc46 = 0, acc56 = 0, acc66 = 0, acc76 = 0;
        register T acc07 = 0, acc17 = 0, acc27 = 0, acc37 = 0, acc47 = 0, acc57 = 0, acc67 = 0, acc77 = 0;

        for (mat_size_t k = 0; k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01 += elements[n_cols * (i + 0) + k] * other(k, j + 1);
            acc02 += elements[n_cols * (i + 0) + k] * other(k, j + 2);
            acc03 += elements[n_cols * (i + 0) + k] * other(k, j + 3);
            acc04 += elements[n_cols * (i + 0) + k] * other(k, j + 4);
            acc05 += elements[n_cols * (i + 0) + k] * other(k, j + 5);
            acc06 += elements[n_cols * (i + 0) + k] * other(k, j + 6);
            acc07 += elements[n_cols * (i + 0) + k] * other(k, j + 7);

            acc10 += elements[n_cols * (i + 1) + k] * other(k, j + 0);
            acc11 += elements[n_cols * (i + 1) + k] * other(k, j + 1);
            acc12 += elements[n_cols * (i + 1) + k] * other(k, j + 2);
            acc13 += elements[n_cols * (i + 1) + k] * other(k, j + 3);
            acc14 += elements[n_cols * (i + 1) + k] * other(k, j + 4);
            acc15 += elements[n_cols * (i + 1) + k] * other(k, j + 5);
            acc16 += elements[n_cols * (i + 1) + k] * other(k, j + 6);
            acc17 += elements[n_cols * (i + 1) + k] * other(k, j + 7);

            acc20 += elements[n_cols * (i + 2) + k] * other(k, j + 0);
            acc21 += elements[n_cols * (i + 2) + k] * other(k, j + 1);
            acc22 += elements[n_cols * (i + 2) + k] * other(k, j + 2);
            acc23 += elements[n_cols * (i + 2) + k] * other(k, j + 3);
            acc24 += elements[n_cols * (i + 2) + k] * other(k, j + 4);
            acc25 += elements[n_cols * (i + 2) + k] * other(k, j + 5);
            acc26 += elements[n_cols * (i + 2) + k] * other(k, j + 6);
            acc27 += elements[n_cols * (i + 2) + k] * other(k, j + 7);

            acc30 += elements[n_cols * (i + 3) + k] * other(k, j + 0);
            acc31 += elements[n_cols * (i + 3) + k] * other(k, j + 1);
            acc32 += elements[n_cols * (i + 3) + k] * other(k, j + 2);
            acc33 += elements[n_cols * (i + 3) + k] * other(k, j + 3);
            acc34 += elements[n_cols * (i + 3) + k] * other(k, j + 4);
            acc35 += elements[n_cols * (i + 3) + k] * other(k, j + 5);
            acc36 += elements[n_cols * (i + 3) + k] * other(k, j + 6);
            acc37 += elements[n_cols * (i + 3) + k] * other(k, j + 7);

            acc40 += elements[n_cols * (i + 4) + k] * other(k, j + 0);
            acc41 += elements[n_cols * (i + 4) + k] * other(k, j + 1);
            acc42 += elements[n_cols * (i + 4) + k] * other(k, j + 2);
            acc43 += elements[n_cols * (i + 4) + k] * other(k, j + 3);
            acc44 += elements[n_cols * (i + 4) + k] * other(k, j + 4);
            acc45 += elements[n_cols * (i + 4) + k] * other(k, j + 5);
            acc46 += elements[n_cols * (i + 4) + k] * other(k, j + 6);
            acc47 += elements[n_cols * (i + 4) + k] * other(k, j + 7);

            acc50 += elements[n_cols * (i + 5) + k] * other(k, j + 0);
            acc51 += elements[n_cols * (i + 5) + k] * other(k, j + 1);
            acc52 += elements[n_cols * (i + 5) + k] * other(k, j + 2);
            acc53 += elements[n_cols * (i + 5) + k] * other(k, j + 3);
            acc54 += elements[n_cols * (i + 5) + k] * other(k, j + 4);
            acc55 += elements[n_cols * (i + 5) + k] * other(k, j + 5);
            acc56 += elements[n_cols * (i + 5) + k] * other(k, j + 6);
            acc57 += elements[n_cols * (i + 5) + k] * other(k, j + 7);

            acc60 += elements[n_cols * (i + 6) + k] * other(k, j + 0);
            acc61 += elements[n_cols * (i + 6) + k] * other(k, j + 1);
            acc62 += elements[n_cols * (i + 6) + k] * other(k, j + 2);
            acc63 += elements[n_cols * (i + 6) + k] * other(k, j + 3);
            acc64 += elements[n_cols * (i + 6) + k] * other(k, j + 4);
            acc65 += elements[n_cols * (i + 6) + k] * other(k, j + 5);
            acc66 += elements[n_cols * (i + 6) + k] * other(k, j + 6);
            acc67 += elements[n_cols * (i + 6) + k] * other(k, j + 7);

            acc70 += elements[n_cols * (i + 7) + k] * other(k, j + 0);
            acc71 += elements[n_cols * (i + 7) + k] * other(k, j + 1);
            acc72 += elements[n_cols * (i + 7) + k] * other(k, j + 2);
            acc73 += elements[n_cols * (i + 7) + k] * other(k, j + 3);
            acc74 += elements[n_cols * (i + 7) + k] * other(k, j + 4);
            acc75 += elements[n_cols * (i + 7) + k] * other(k, j + 5);
            acc76 += elements[n_cols * (i + 7) + k] * other(k, j + 6);
            acc77 += elements[n_cols * (i + 7) + k] * other(k, j + 7);
        }  // k
        res(i + 0, j + 0) = acc00;
        res(i + 0, j + 1) = acc01;
        res(i + 0, j + 2) = acc02;
        res(i + 0, j + 3) = acc03;
        res(i + 0, j + 4) = acc04;
        res(i + 0, j + 5) = acc05;
        res(i + 0, j + 6) = acc06;
        res(i + 0, j + 7) = acc07;

        res(i + 1, j + 0) = acc10;
        res(i + 1, j + 1) = acc11;
        res(i + 1, j + 2) = acc12;
        res(i + 1, j + 3) = acc13;
        res(i + 1, j + 4) = acc14;
        res(i + 1, j + 5) = acc15;
        res(i + 1, j + 6) = acc16;
        res(i + 1, j + 7) = acc17;

        res(i + 2, j + 0) = acc20;
        res(i + 2, j + 1) = acc21;
        res(i + 2, j + 2) = acc22;
        res(i + 2, j + 3) = acc23;
        res(i + 2, j + 4) = acc24;
        res(i + 2, j + 5) = acc25;
        res(i + 2, j + 6) = acc26;
        res(i + 2, j + 7) = acc27;

        res(i + 3, j + 0) = acc30;
        res(i + 3, j + 1) = acc31;
        res(i + 3, j + 2) = acc32;
        res(i + 3, j + 3) = acc33;
        res(i + 3, j + 4) = acc34;
        res(i + 3, j + 5) = acc35;
        res(i + 3, j + 6) = acc36;
        res(i + 3, j + 7) = acc37;

        res(i + 4, j + 0) = acc40;
        res(i + 4, j + 1) = acc41;
        res(i + 4, j + 2) = acc42;
        res(i + 4, j + 3) = acc43;
        res(i + 4, j + 4) = acc44;
        res(i + 4, j + 5) = acc45;
        res(i + 4, j + 6) = acc46;
        res(i + 4, j + 7) = acc47;

        res(i + 5, j + 0) = acc50;
        res(i + 5, j + 1) = acc51;
        res(i + 5, j + 2) = acc52;
        res(i + 5, j + 3) = acc53;
        res(i + 5, j + 4) = acc54;
        res(i + 5, j + 5) = acc55;
        res(i + 5, j + 6) = acc56;
        res(i + 5, j + 7) = acc57;

        res(i + 6, j + 0) = acc60;
        res(i + 6, j + 1) = acc61;
        res(i + 6, j + 2) = acc62;
        res(i + 6, j + 3) = acc63;
        res(i + 6, j + 4) = acc64;
        res(i + 6, j + 5) = acc65;
        res(i + 6, j + 6) = acc66;
        res(i + 6, j + 7) = acc67;

        res(i + 7, j + 0) = acc70;
        res(i + 7, j + 1) = acc71;
        res(i + 7, j + 2) = acc72;
        res(i + 7, j + 3) = acc73;
        res(i + 7, j + 4) = acc74;
        res(i + 7, j + 5) = acc75;
        res(i + 7, j + 6) = acc76;
        res(i + 7, j + 7) = acc77;
    }

    inline void blockDotNx1(mat_size_t i, mat_size_t j, Matrix<T> &other, Matrix<T> &res) {
        register T acc00a = 0, acc10a = 0, acc20a = 0, acc30a = 0;
        register T acc40a = 0, acc50a = 0, acc60a = 0, acc70a = 0;
        register T acc00b = 0, acc10b = 0, acc20b = 0, acc30b = 0;
        register T acc40b = 0, acc50b = 0, acc60b = 0, acc70b = 0;
        mat_size_t k;
        for (k = 0; k + 2 < this->n_cols; k += 2) {
            acc00a += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc10a += elements[n_cols * (i + 1) + k] * other(k, j + 0);
            acc20a += elements[n_cols * (i + 2) + k] * other(k, j + 0);
            acc30a += elements[n_cols * (i + 3) + k] * other(k, j + 0);
            acc40a += elements[n_cols * (i + 4) + k] * other(k, j + 0);
            acc50a += elements[n_cols * (i + 5) + k] * other(k, j + 0);
            acc60a += elements[n_cols * (i + 6) + k] * other(k, j + 0);
            acc70a += elements[n_cols * (i + 7) + k] * other(k, j + 0);

            acc00b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 0);
            acc10b += elements[n_cols * (i + 1) + k + 1] * other(k + 1, j + 0);
            acc20b += elements[n_cols * (i + 2) + k + 1] * other(k + 1, j + 0);
            acc30b += elements[n_cols * (i + 3) + k + 1] * other(k + 1, j + 0);
            acc40b += elements[n_cols * (i + 4) + k + 1] * other(k + 1, j + 0);
            acc50b += elements[n_cols * (i + 5) + k + 1] * other(k + 1, j + 0);
            acc60b += elements[n_cols * (i + 6) + k + 1] * other(k + 1, j + 0);
            acc70b += elements[n_cols * (i + 7) + k + 1] * other(k + 1, j + 0);
        }  // k
        for (; k < this->n_cols; ++k) {
            acc00a += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc10a += elements[n_cols * (i + 1) + k] * other(k, j + 0);
            acc20a += elements[n_cols * (i + 2) + k] * other(k, j + 0);
            acc30a += elements[n_cols * (i + 3) + k] * other(k, j + 0);
            acc40a += elements[n_cols * (i + 4) + k] * other(k, j + 0);
            acc50a += elements[n_cols * (i + 5) + k] * other(k, j + 0);
            acc60a += elements[n_cols * (i + 6) + k] * other(k, j + 0);
            acc70a += elements[n_cols * (i + 7) + k] * other(k, j + 0);
        }

        res(i + 0, j + 0) = acc00a + acc00b;
        res(i + 1, j + 0) = acc10a + acc10b;
        res(i + 2, j + 0) = acc20a + acc20b;
        res(i + 3, j + 0) = acc30a + acc30b;
        res(i + 4, j + 0) = acc40a + acc40b;
        res(i + 5, j + 0) = acc50a + acc50b;
        res(i + 6, j + 0) = acc60a + acc60b;
        res(i + 7, j + 0) = acc70a + acc70b;
    }

    inline void blockDot1xN(mat_size_t i, mat_size_t j, Matrix<T> &other, Matrix<T> &res) {
        register T acc00a = 0, acc01a = 0, acc02a = 0, acc03a = 0;
        register T acc04a = 0, acc05a = 0, acc06a = 0, acc07a = 0;
        register T acc00b = 0, acc01b = 0, acc02b = 0, acc03b = 0;
        register T acc04b = 0, acc05b = 0, acc06b = 0, acc07b = 0;
        mat_size_t k;
        for (k = 0; k + 2 < this->n_cols; k += 2) {
            acc00a += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01a += elements[n_cols * (i + 0) + k] * other(k, j + 1);
            acc02a += elements[n_cols * (i + 0) + k] * other(k, j + 2);
            acc03a += elements[n_cols * (i + 0) + k] * other(k, j + 3);
            acc04a += elements[n_cols * (i + 0) + k] * other(k, j + 4);
            acc05a += elements[n_cols * (i + 0) + k] * other(k, j + 5);
            acc06a += elements[n_cols * (i + 0) + k] * other(k, j + 6);
            acc07a += elements[n_cols * (i + 0) + k] * other(k, j + 7);

            acc00b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 0);
            acc01b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 1);
            acc02b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 2);
            acc03b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 3);
            acc04b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 4);
            acc05b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 5);
            acc06b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 6);
            acc07b += elements[n_cols * (i + 0) + k + 1] * other(k + 1, j + 7);
        }  // k
        for (; k < this->n_cols; ++k) {
            acc00a += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01a += elements[n_cols * (i + 0) + k] * other(k, j + 1);
            acc02a += elements[n_cols * (i + 0) + k] * other(k, j + 2);
            acc03a += elements[n_cols * (i + 0) + k] * other(k, j + 3);
            acc04a += elements[n_cols * (i + 0) + k] * other(k, j + 4);
            acc05a += elements[n_cols * (i + 0) + k] * other(k, j + 5);
            acc06a += elements[n_cols * (i + 0) + k] * other(k, j + 6);
            acc07a += elements[n_cols * (i + 0) + k] * other(k, j + 7);
        }  // k
        res(i + 0, j + 0) = acc00a + acc00b;
        res(i + 0, j + 1) = acc01a + acc01b;
        res(i + 0, j + 2) = acc02a + acc02b;
        res(i + 0, j + 3) = acc03a + acc03b;
        res(i + 0, j + 4) = acc04a + acc04b;
        res(i + 0, j + 5) = acc05a + acc05b;
        res(i + 0, j + 6) = acc06a + acc06b;
        res(i + 0, j + 7) = acc07a + acc07b;
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
