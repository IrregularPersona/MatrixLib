#ifndef MATRIX_H
#define MATRIX_H

#include <array>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <sstream>

#define INT_DATATYPE int64_t
// typedef __int128_t INT_DATATYPE;

template <size_t Rows, size_t Cols>
class Matrix {
private:
    std::array<INT_DATATYPE, Rows * Cols> m_Data;

public:
    Matrix();

    // void print128 (INT_DATATYPE val);

    void set(size_t i, size_t j);
    const INT_DATATYPE& get(size_t i, size_t j) const;

    void identity();
    void debug_print() const;

    Matrix& init(INT_DATATYPE value);

    Matrix<Cols, Rows> transpose() const;
    Matrix<Rows, Cols> inverse() const;
    
    INT_DATATYPE determinant() const;
    Matrix<Rows - 1, Cols - 1> minor_matrix(size_t row, size_t col) const;


    Matrix operator+(const Matrix& other) const;

    Matrix operator-(const Matrix& other) const;    
    Matrix operator-() const;

    template <size_t OtherCols>
    Matrix<Rows, OtherCols> operator*(const Matrix<Cols, OtherCols>& other) const;

    Matrix operator*(INT_DATATYPE scalar) const;
    
    void swapRows(size_t row1, size_t row2);

    size_t size() const;
};

template <size_t Rows, size_t Cols>
Matrix<Rows, Cols> operator*(INT_DATATYPE scalar, const Matrix<Rows, Cols>& matrix);

constexpr bool isPerfectSquare(size_t n);

template<size_t N>
INT_DATATYPE determinantLU(const Matrix<N, N>& mat);

#include "matrix.tpp" // Include the template implementation
#endif // MATRIX_H
