#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

template <class T>
class Matrix {
public:
    explicit Matrix(std::vector<std::vector<T> > _elements) : elements(_elements) {
        setShape();
    };

    Matrix(size_t n_rows, size_t n_cols) :
            Matrix(std::vector<std::vector<T> >(n_rows, std::vector<T>(n_cols))) {}

    Matrix() : Matrix(std::vector<std::vector<T> >()) {}

    bool empty() const {
        return elements.empty();
    }

    Matrix transpose() const {
        Matrix<T> mT(_shape[1], _shape[0]);
        for (size_t i = 0; i < _shape[1]; ++i) {
            for (size_t j = 0; j < _shape[0]; ++j)
                mT[i][j] = elements[j][i];
        }

        return mT;
    }

    std::string to_string() const {
        std::stringstream ss;
        for (std::vector<T> row : elements) {
            for (T elem : row)
                ss << elem << '\t';
            ss << '\n';
        }
        return ss.str();
    }
    size_t shape(size_t i) const {
        return _shape[i];
    }

    std::vector<T>& operator[](size_t i) {
        return elements[i];
    }

    Matrix<T> operator*(Matrix<T>& m) const {
        assert(!(empty() || m.empty()));
        assert(_shape[1] == m.shape(0));
        Matrix<T> res(_shape[0], m.shape(1));
        for (size_t i = 0; i < _shape[0]; ++i)
            for (size_t k = 0; k < m.shape(1); ++k)
                for (size_t j = 0; j < _shape[1]; ++j)
                    res[i][k] += elements[i][j] * m[j][k];
        return res;
    }

private:
    std::vector<std::vector<T> > elements;
    std::vector<size_t> _shape;

    void setShape() {
        _shape = std::vector<size_t>(2);
        if (!empty()) {
            _shape[0] = elements.size();
            _shape[1] = elements[0].size();
        }
    }
};

template<class T>
static std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {
    return os << m.to_string();
}

#endif //MATRIX_MATRIX_H
