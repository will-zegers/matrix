#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <vector>
#include <cstdint>
#include <sstream>
#include <tuple>
#include <iostream>

#define XPOSE_STEP 2
#define MATMUL_STEP 8
#define I_BLOCK_SZ 64
#define K_BLOCK_SZ 32

typedef uint32_t mat_size_t;
typedef std::tuple<mat_size_t, mat_size_t> shape_t;

template <typename T>
class Matrix {
public:
    /**
     * Instatiates an empty matrix of size (0 x 0)
     */
    Matrix() : n_rows(0), n_cols(0), elements(std::vector<T> ()) {}
    /**
     * @param shape Tuple containg the number of rows in the first element and number of columns in the second
     */
    Matrix(shape_t shape) :
            n_rows(std::get<0>(shape)),
            n_cols(std::get<1>(shape)),
            elements(std::vector<T> (n_rows * n_cols)) {}
    /**
     * @param shape Tuple containg the number of rows in the first element and number of columns in the second
     *  @param elements A standard vector containg the m x n elements of the matrix.
     */
    Matrix(shape_t shape, std::vector<T> elements) :
            n_rows(std::get<0>(shape)),
            n_cols(std::get<1>(shape)),
            elements(elements) {}

    /**
     * @param i Selected row
     * @param j Selected column
     * @return The element at index [i, j]
     */
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

    /**
     * Return a copy of this matrix instance, transposed. Uses loop unrolling to gain a little bit of a performance
     * boost
     * @return A new Matrix instance.
     */
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

    /**
     * Optimized matrix multiplication mostly drawn from Wikipedia: https://goo.gl/JRWcsB, and some conversations on
     * StackOverflow. Uses loop tiling and unrolling to perform matrix multiplication.
     *
     * @param other Another matrix instance.
     * @return A Matrix instance resulting from the multiplication of this and other.
     */
    virtual Matrix<T> operator*(Matrix<T>& other) {
        if (this->empty() || other.empty())
            throw empty_matrix();
        if (this->n_cols != other.shape(0)) {
            throw size_mismatch();
        }

        Matrix<T> res = Matrix<T>(std::make_pair(this->n_rows, other.shape(1)));
        mat_size_t ii;
        for (ii = 0; ii + I_BLOCK_SZ < this->n_rows; ii += I_BLOCK_SZ) {
            mat_size_t kk;
            for (kk = 0; kk + K_BLOCK_SZ < this->n_cols; kk += K_BLOCK_SZ) {
                mat_size_t j;
                for (j = 0; j + MATMUL_STEP < other.shape(1); j += MATMUL_STEP) {
                    mat_size_t i;
                    for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                        this->blockDotNxN(i, j, kk, other, res);
                    }  // i
                    for (; i < ii + I_BLOCK_SZ; ++i) {  // Clean up the last few rows that didn't align with MATMUL_STEP
                        this->blockDot1xN(i, j, kk, other, res);
                    }  // i
                }  // j
                for (; j < other.shape(1); ++j) {  // Clean up the last few columns that didn't align with MATMUL_STEP
                    mat_size_t i;
                    for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                        this->blockDotNx1(i, j, kk, other, res);
                    }  // i
                    for (; i < ii + I_BLOCK_SZ; ++i) {  // Clean up the last few rows that didn't align with MATMUL_STEP
                        this->dot(i, j, kk, other, res);
                    }  // i
                }  // j
            }  // kk
            mat_size_t j;  // Clean up columns and row that didn't align with K_BLOCK_SZ
            for (j = 0; j + MATMUL_STEP < other.shape(1); j += MATMUL_STEP) {
                mat_size_t i;
                for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                    this->blockDotNxN(i, j, kk, other, res);
                }  // i
                for (; i < ii + I_BLOCK_SZ; ++i) {  // Clean up the last few rows that didn't align with MATMUL_STEP
                    this->blockDot1xN(i, j, kk, other, res);
                }  // i
            }  // j
            for (; j < other.shape(1); ++j) {  // Clean up the last few columns that didn't align with MATMUL_STEP
                mat_size_t i;
                for (i = ii; i < ii + I_BLOCK_SZ; i += MATMUL_STEP) {
                    this->blockDotNx1(i, j, kk, other, res);
                }  // i
                for (; i < ii + I_BLOCK_SZ; ++i) {  // Clean up the last few rows that didn't align with MATMUL_STEP
                    this->dot(i, j, kk, other, res);
                }  // i
            }  // j
        }  // ii
        mat_size_t kk;  // Clean up columns and row that didn't align with I_BLOCK_SZ
        for (kk = 0; kk + K_BLOCK_SZ < this->n_cols; kk += K_BLOCK_SZ) {
            mat_size_t j;
            for (j = 0; j + MATMUL_STEP < other.shape(1); j += MATMUL_STEP) {
                mat_size_t i;
                for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                    this->blockDotNxN(i, j, kk, other, res);
                }  // i
                for (; i < this->n_rows; ++i) {  // Clean up the last few rows that didn't align with MATMUL_STEP
                    this->blockDot1xN(i, j, kk, other, res);
                }  // i
            }  // j
            for (; j < other.shape(1); ++j) {  // Clean up the last few columns that didn't align with MATMUL_STEP
                mat_size_t i;
                for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                    this->blockDotNx1(i, j, kk, other, res);
                }  // i
                for (; i < this->n_rows; ++i) {  // Clean up the last few rows that didn't align with MATMUL_STEP
                    this->dot(i, j, kk, other, res);
                }  // i
            }  // j
        }  // kk
        mat_size_t j;  // Clean up columns and row that didn't align with K_BLOCK_SZ
        for (j = 0; j + MATMUL_STEP < other.shape(1); j += MATMUL_STEP) {
            mat_size_t i;
            for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                this->blockDotNxN(i, j, kk, other, res);
            }  // i
            for (; i < this->n_rows; ++i) {  // Clean up the last few rows
                this->blockDot1xN(i, j, kk, other, res);
            }  // i
        }  // j
        for (; j < other.shape(1); ++j) {  // Clean up the remaining last few columns
            mat_size_t i;
            for (i = ii; i + MATMUL_STEP < this->n_rows; i += MATMUL_STEP) {
                this->blockDotNx1(i, j, kk,  other, res);
            }  // i
            for (; i < this->n_rows; ++i) {  // Clean up the last few rows
                this->dot(i, j, kk, other, res);
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

    /**
     * @param n 0 for number of rows, 1 for number of columns
     * @return The requested dimension
     */
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

    /**
     * Thrown when the matrix has dimensions 0x0 (empty)
     */
    struct empty_matrix : public std::exception {
        const char* what() const throw() final {
            return "Cannot perform operation on an empty matrix";
        }
    };

    /**
     * Thrown when trying to multiply (m x n) and (n' x m') matrix, where n != n'
     */
    struct size_mismatch : public std::exception {
        const char* what() const throw() final {
            return "Matrix dimensions do not match";
        }
    };

    /**
     * Thrown when requesting an invalid shape parameter (i.e. anything other than 0 or 1)
     */
    struct bad_shape : public std::exception {
        const char* what() const throw() final {
            return "Requested dimension does not exist";
        }
    };

protected:
    mat_size_t n_rows, n_cols;
    std::vector<T> elements;

    /**
     * Computes the dot product over NxN blocks of elements of two matrices
     *
     * @param i Current index of i in the matrix multiplication operation
     * @param j Current index of j in the matrix multiplication operation
     * @param kk Current index of kk in the matrix multiplication operation
     * @param other An instance of the multiplicand(?) Matrix
     * @param res The result which into which the dot product will be stored
     */
    inline void blockDotNxN(mat_size_t i, mat_size_t j, mat_size_t kk, Matrix<T> &other, Matrix<T> &res) {
        register T acc00, acc10, acc20, acc30, acc40, acc50, acc60, acc70;
        register T acc01, acc11, acc21, acc31, acc41, acc51, acc61, acc71;
        register T acc02, acc12, acc22, acc32, acc42, acc52, acc62, acc72;
        register T acc03, acc13, acc23, acc33, acc43, acc53, acc63, acc73;
        register T acc04, acc14, acc24, acc34, acc44, acc54, acc64, acc74;
        register T acc05, acc15, acc25, acc35, acc45, acc55, acc65, acc75;
        register T acc06, acc16, acc26, acc36, acc46, acc56, acc66, acc76;
        register T acc07, acc17, acc27, acc37, acc47, acc57, acc67, acc77;

        if (kk == 0) {
            acc00 = 0, acc10 = 0, acc20 = 0, acc30 = 0, acc40 = 0, acc50 = 0, acc60 = 0, acc70 = 0;
            acc01 = 0, acc11 = 0, acc21 = 0, acc31 = 0, acc41 = 0, acc51 = 0, acc61 = 0, acc71 = 0;
            acc02 = 0, acc12 = 0, acc22 = 0, acc32 = 0, acc42 = 0, acc52 = 0, acc62 = 0, acc72 = 0;
            acc03 = 0, acc13 = 0, acc23 = 0, acc33 = 0, acc43 = 0, acc53 = 0, acc63 = 0, acc73 = 0;
            acc04 = 0, acc14 = 0, acc24 = 0, acc34 = 0, acc44 = 0, acc54 = 0, acc64 = 0, acc74 = 0;
            acc05 = 0, acc15 = 0, acc25 = 0, acc35 = 0, acc45 = 0, acc55 = 0, acc65 = 0, acc75 = 0;
            acc06 = 0, acc16 = 0, acc26 = 0, acc36 = 0, acc46 = 0, acc56 = 0, acc66 = 0, acc76 = 0;
            acc07 = 0, acc17 = 0, acc27 = 0, acc37 = 0, acc47 = 0, acc57 = 0, acc67 = 0, acc77 = 0;
        } else {
            acc00 = res(i + 0, j + 0); acc04 = res(i + 0, j + 4);
            acc01 = res(i + 0, j + 1); acc05 = res(i + 0, j + 5);
            acc02 = res(i + 0, j + 2); acc06 = res(i + 0, j + 6);
            acc03 = res(i + 0, j + 3); acc07 = res(i + 0, j + 7);

            acc10 = res(i + 1, j + 0); acc14 = res(i + 1, j + 4);
            acc11 = res(i + 1, j + 1); acc15 = res(i + 1, j + 5);
            acc12 = res(i + 1, j + 2); acc16 = res(i + 1, j + 6);
            acc13 = res(i + 1, j + 3); acc17 = res(i + 1, j + 7);

            acc20 = res(i + 2, j + 0); acc24 = res(i + 2, j + 4);
            acc21 = res(i + 2, j + 1); acc25 = res(i + 2, j + 5);
            acc22 = res(i + 2, j + 2); acc26 = res(i + 2, j + 6);
            acc23 = res(i + 2, j + 3); acc27 = res(i + 2, j + 7);

            acc30 = res(i + 3, j + 0); acc34 = res(i + 3, j + 4);
            acc31 = res(i + 3, j + 1); acc35 = res(i + 3, j + 5);
            acc32 = res(i + 3, j + 2); acc36 = res(i + 3, j + 6);
            acc33 = res(i + 3, j + 3); acc37 = res(i + 3, j + 7);

            acc40 = res(i + 4, j + 0); acc44 = res(i + 4, j + 4);
            acc41 = res(i + 4, j + 1); acc45 = res(i + 4, j + 5);
            acc42 = res(i + 4, j + 2); acc46 = res(i + 4, j + 6);
            acc43 = res(i + 4, j + 3); acc47 = res(i + 4, j + 7);

            acc50 = res(i + 5, j + 0); acc54 = res(i + 5, j + 4);
            acc51 = res(i + 5, j + 1); acc55 = res(i + 5, j + 5);
            acc52 = res(i + 5, j + 2); acc56 = res(i + 5, j + 6);
            acc53 = res(i + 5, j + 3); acc57 = res(i + 5, j + 7);

            acc60 = res(i + 6, j + 0); acc64 = res(i + 6, j + 4);
            acc61 = res(i + 6, j + 1); acc65 = res(i + 6, j + 5);
            acc62 = res(i + 6, j + 2); acc66 = res(i + 6, j + 6);
            acc63 = res(i + 6, j + 3); acc67 = res(i + 6, j + 7);

            acc70 = res(i + 7, j + 0); acc74 = res(i + 7, j + 4);
            acc71 = res(i + 7, j + 1); acc75 = res(i + 7, j + 5);
            acc72 = res(i + 7, j + 2); acc76 = res(i + 7, j + 6);
            acc73 = res(i + 7, j + 3); acc77 = res(i + 7, j + 7);
        }

        for (mat_size_t k = kk; k < kk + K_BLOCK_SZ && k < this->n_cols; ++k) {
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

    /**
     * Computes the dot product over Nx1 blocks of elements of two matrices. Mainly used to clean up the last few
     * columns along a block of multilpications that don't align with the block.
     *
     * @param i Current index of i in the matrix multiplcation operation
     * @param j Current index of j in the matrix multiplcation operation
     * @param kk Current index of kk in the matrix multiplcation operation
     * @param other An instance of the multiplicand(?) Matrix
     * @param res The result which into which the dot product will be stored
     */
    inline void blockDotNx1(mat_size_t i, mat_size_t j, mat_size_t kk, Matrix<T> &other, Matrix<T> &res) {
        register T acc00, acc10, acc20, acc30, acc40, acc50, acc60, acc70;

        if (kk == 0) {
            acc00 = 0, acc10 = 0, acc20 = 0, acc30 = 0, acc40 = 0, acc50 = 0, acc60 = 0, acc70 = 0;
        } else {
            acc00 = res(i + 0, j + 0); acc40 = res(i + 4, j + 0);
            acc10 = res(i + 1, j + 0); acc50 = res(i + 5, j + 0);
            acc20 = res(i + 2, j + 0); acc60 = res(i + 6, j + 0);
            acc30 = res(i + 3, j + 0); acc70 = res(i + 7, j + 0);
        }

        for (mat_size_t k = kk; k < kk + K_BLOCK_SZ && k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc10 += elements[n_cols * (i + 1) + k] * other(k, j + 0);
            acc20 += elements[n_cols * (i + 2) + k] * other(k, j + 0);
            acc30 += elements[n_cols * (i + 3) + k] * other(k, j + 0);
            acc40 += elements[n_cols * (i + 4) + k] * other(k, j + 0);
            acc50 += elements[n_cols * (i + 5) + k] * other(k, j + 0);
            acc60 += elements[n_cols * (i + 6) + k] * other(k, j + 0);
            acc70 += elements[n_cols * (i + 7) + k] * other(k, j + 0);
        }  // k

        res(i + 0, j + 0) = acc00;
        res(i + 1, j + 0) = acc10;
        res(i + 2, j + 0) = acc20;
        res(i + 3, j + 0) = acc30;
        res(i + 4, j + 0) = acc40;
        res(i + 5, j + 0) = acc50;
        res(i + 6, j + 0) = acc60;
        res(i + 7, j + 0) = acc70;
    }

    /**
     * Computes the dot product over Nx1 blocks of elements of two matrices. Mainly used to clean up the last few
     * rows along a block of multiplications that don't align with the block.
     *
     * @param i Current index of i in the matrix multiplication operation
     * @param j Current index of j in the matrix multiplication operation
     * @param kk Current index of kk in the matrix multiplication operation
     * @param other An instance of the multiplicand(?) Matrix
     * @param res The result which into which the dot product will be stored
     */
    inline void blockDot1xN(mat_size_t i, mat_size_t j, mat_size_t kk, Matrix<T> &other, Matrix<T> &res) {
        register T acc00, acc01, acc02, acc03, acc04, acc05, acc06, acc07;

        if (kk == 0) {
            acc00 = 0;  acc01 = 0; acc02 = 0; acc03 = 0; acc04 = 0; acc05 = 0; acc06 = 0; acc07 = 0;
        } else {
            acc00 = res(i + 0, j + 0); acc04 = res(i + 0, j + 4);
            acc01 = res(i + 0, j + 1); acc05 = res(i + 0, j + 5);
            acc02 = res(i + 0, j + 2); acc06 = res(i + 0, j + 6);
            acc03 = res(i + 0, j + 3); acc07 = res(i + 0, j + 7);
        }

        for (mat_size_t k = kk; k < kk + K_BLOCK_SZ && k < this->n_cols; ++k) {
            acc00 += elements[n_cols * (i + 0) + k] * other(k, j + 0);
            acc01 += elements[n_cols * (i + 0) + k] * other(k, j + 1);
            acc02 += elements[n_cols * (i + 0) + k] * other(k, j + 2);
            acc03 += elements[n_cols * (i + 0) + k] * other(k, j + 3);
            acc04 += elements[n_cols * (i + 0) + k] * other(k, j + 4);
            acc05 += elements[n_cols * (i + 0) + k] * other(k, j + 5);
            acc06 += elements[n_cols * (i + 0) + k] * other(k, j + 6);
            acc07 += elements[n_cols * (i + 0) + k] * other(k, j + 7);
        }  // k

        res(i + 0, j + 0) = acc00;
        res(i + 0, j + 1) = acc01;
        res(i + 0, j + 2) = acc02;
        res(i + 0, j + 3) = acc03;
        res(i + 0, j + 4) = acc04;
        res(i + 0, j + 5) = acc05;
        res(i + 0, j + 6) = acc06;
        res(i + 0, j + 7) = acc07;
    }

    /**
     * Computes just the ordinary dot product. Used to clean up the last few elements that well outside the row and cell
     * blocks.
     *
     * @param i Current index of i in the matrix multiplication operation
     * @param j Current index of j in the matrix multiplication operation
     * @param kk Current index of kk in the matrix multiplication operation
     * @param other An instance of the multiplicand(?) Matrix
     * @param res The result which into which the dot product will be stored
     */
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
