#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <iostream>
#include <exception>

typedef uint32_t mat_size_t;

template <typename T>
class Matrix {
public:
    typedef std::vector<std::vector<T> > vector2D;

    Matrix() : Matrix(vector2D()) {}

    Matrix(mat_size_t n_rows, mat_size_t n_cols) :
            Matrix(vector2D(n_rows, std::vector<T>(n_cols))) {}

    explicit Matrix(vector2D _elements) : elements(_elements) { setShape(); };

    typename
    vector2D::iterator begin() { return elements.begin(); }

    typename
    vector2D::iterator end() { return elements.end(); }

    bool empty() const { return elements.empty(); }

    mat_size_t shape(mat_size_t i) const { return _shape[i]; }

    std::vector<T> at(mat_size_t i) const { return elements[i]; }

    T at(mat_size_t i, mat_size_t j) const { return elements[i][j]; }

    std::string to_string() const {
        std::stringstream ss;
        for (std::vector<T> row : elements) {
            for (T elem : row)
                ss << elem << '\t';
            ss << '\n';
        }
        return ss.str();
    }

    std::vector<T>& operator[](mat_size_t i) { return elements[i]; }

    bool operator==(const Matrix& m) const {
        if (_shape[0] != m.shape(0) || _shape[1] != m.shape(1))
            return false;

        for (mat_size_t i = 0; i < _shape[0]; ++i)
            if (elements.at(i) != m.at(i))
                return false;
        return true;
    }

    virtual Matrix<T> operator*(Matrix<T>& m) const {
        if (empty() || m.empty())
            throw empty_matrix();
        if (_shape[1] != m.shape(0))
            throw size_mismatch();

        Matrix<T> res(_shape[0], m.shape(1));
        for (mat_size_t i = 0; i < _shape[0]; ++i)
            for (mat_size_t k = 0; k < m.shape(1); ++k)
                for (mat_size_t j = 0; j < _shape[1]; ++j)
                    res[i][k] += elements[i][j] * m[j][k];
        return res;
    }

    virtual Matrix transpose() {
        if (empty())
            throw empty_matrix();

        Matrix<T> mT(_shape[1], _shape[0]);
        for (mat_size_t i = 0; i < _shape[1]; ++i)
            for (mat_size_t j = 0; j < _shape[0]; ++j)
                mT[i][j] = elements[j][i];

        return mT;
    }

    struct empty_matrix : public std::exception {
        const char* what () const throw() final {
            return "Cannot perform operation on empty matrices";
        }
    };

    struct size_mismatch : public std::exception {
        const char* what () const throw() final {
            return "Incompatible matrix dimensions";
        }
    };

protected:
    vector2D elements;
    std::vector<mat_size_t> _shape;

    void setShape() {
        _shape = std::vector<mat_size_t>(2);
        if (!empty()) {
            _shape[0] = elements.size();
            _shape[1] = elements[0].size();
        }
    }
};

template <typename T>
static std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {
    return os << m.to_string();
}

#endif //MATRIX_MATRIX_H
