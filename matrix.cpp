#include <iostream>
#include "matrix.h"
#include <sstream>
#include <cassert>

template<class T>
Matrix<T>::Matrix(std::vector<std::vector<T> > _elements) : elements(_elements) {
    Matrix<T>::setShape();
}

template<class T>
Matrix<T>::Matrix(size_t n_rows, size_t n_cols) :
        Matrix<T>::Matrix(std::vector<std::vector<T> >(
                n_rows,
                std::vector<T>(n_cols))) {}

template<class T>
Matrix<T>::Matrix() : elements(std::vector<std::vector<T> >()) {}

template<class T>
bool Matrix<T>::isEmpty() const {
    return elements.begin() == elements.end();
}

template<class T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> mT(_shape[1], _shape[0]);
    for (size_t i = 0; i < _shape[1]; ++i) {
        for (size_t j = 0; j < _shape[0]; ++j)
            mT[i][j] = elements[j][i];
    }

    return mT;
}

template<class T>
std::string Matrix<T>::toString() const {
    std::stringstream ss;
    for (std::vector<T> row : elements) {
        for (T elem : row)
            ss << elem << "\t";
        ss << "\n";
    }
    return ss.str();
}

template<class T>
size_t Matrix<T>::shape(size_t i) const {
    return _shape[i];
}

template<class T>
std::vector<T>& Matrix<T>::operator[](size_t i) {
    return elements[i];
}

template<class T>
T Matrix<T>::at(size_t i, size_t j) const {
    return elements[i][j];
}

template<class T>
void Matrix<T>::setShape() {
    _shape = std::vector<size_t>(2);
    _shape[0] = elements.size();
    _shape[1] = elements[0].size();
}

template<class T>
static std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {
    return os << m.toString();
}

template<class T>
static Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2) {
    assert(!(m1.isEmpty() || m1.isEmpty()));
    assert(m1.shape(1) == m2.shape(0));
    Matrix<T> res(m1.shape(0), m2.shape(1));
    for (size_t i = 0; i < m1.shape(0); ++i)
        for (size_t k = 0; k < m2.shape(1); ++k)
            for (size_t j = 0; j < m1.shape(1); ++j)
                res[i][k] += m1.at(i, j) * m2.at(j, k);
    return res;
}