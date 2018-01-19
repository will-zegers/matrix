#include <cstdlib>
#include <vector>
#include <string>

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

template<class T>
class Matrix {
public:
    explicit Matrix(std::vector<std::vector<T> >);
    explicit Matrix(std::vector<T>);
    Matrix(size_t, size_t);
    Matrix();

    bool isEmpty() const;
    Matrix transpose() const;
    std::string toString() const;
    size_t shape(size_t) const;
    std::vector<T>& operator[](size_t);

private:
    std::vector<std::vector<T> > elements;
    std::vector<size_t> _shape;

    void setShape();
};

template<class T>
static std::ostream& operator<<(std::ostream&, const Matrix<T>&);

template<class T>
static Matrix<T> operator*(Matrix<T>&, Matrix<T>&);

#endif //MATRIX_MATRIX_H
