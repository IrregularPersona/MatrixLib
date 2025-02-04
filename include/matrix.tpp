#ifndef MATRIX_TPP
#define MATRIX_TPP

#include "matrix.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <ratio>
#include <stdexcept>
#include <type_traits>
#define _SIZE_CHECK static_assert(Rows == Cols, "Matrix must be square to calculate determinant!");


// void print128 (INT_DATATYPE val) {
// 	unsigned long long high = val >> 64;
// 	unsigned long long low = val & 0xFFFFFFFFFFFFFFFF;
// 	printf("%llu%llu\n", high, low);
// }

template <size_t Rows, size_t Cols>
Matrix<Rows, Cols>::Matrix() {
    m_Data.fill(0); // Default initialize the matrix with zeros
}

/* 
    FOR DETERMINANT / ANY SQUARE MATRIX CHECKS 
    NOT ENTIRELY NECESSARY IG? 
*/
constexpr bool isPerfectSquare(size_t n) {
    if (n == 0 || n == 1) return true;

    int lastBits = n & 0xF;
    if (lastBits !=  0 && lastBits != 1 && lastBits != 4 && lastBits != 9)
        return false;

    // Bin search
    size_t left = 0, right = n;
    while (left <= right) {
        size_t mid = left + (right - left) / 2;
        size_t square = mid * mid;

        if (square == n) return true;
        if (square < n)
            left = mid + 1;
        else
            right = mid - 1;
    }

    return false;
}

/* SET VALUE */
template <size_t Rows, size_t Cols>
INT_DATATYPE& Matrix<Rows, Cols>::set(size_t i, size_t j) {
    if (i >= Rows || j >= Cols) {
        throw std::out_of_range("Index out of range");
    }
    return m_Data[i * Cols + j];
}

/* GET VALUE */
template <size_t Rows, size_t Cols>
const INT_DATATYPE& Matrix<Rows, Cols>::get(size_t i, size_t j) const {
    if (i >= Rows || j >= Cols) {
        throw std::out_of_range("Index out of range");
    }
    return m_Data[i * Cols + j];
}

/* SET MATRIX TO ITS IDENTITY */
template <size_t Rows, size_t Cols>
void Matrix<Rows, Cols>::identity() {
    m_Data.fill(0);
    for (size_t i = 0; i < std::min(Rows, Cols); ++i) {
        set(i, i) = 1;
    }
}

/* __DEBUG TOOL__ */
template <size_t Rows, size_t Cols>
void Matrix<Rows, Cols>::debug_print() const {
    for (size_t i = 0; i < Rows; ++i) {
        for (size_t j = 0; j < Cols; ++j) {
            std::cout << get(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

/* INITIAL VALUE SET */
template <size_t Rows, size_t Cols>
Matrix<Rows, Cols>& Matrix<Rows, Cols>::init(INT_DATATYPE value) {
    m_Data.fill(value);
    return *this;
}


/* RETURNS A TRANSPOSED MATRIX OBJECT */
template <size_t Rows, size_t Cols>
Matrix<Cols, Rows> Matrix<Rows, Cols>::transpose() const {
    Matrix<Cols, Rows> transposed;

    for (size_t i = 0; i < Rows; ++i) {
        for (size_t j = 0; j < Cols; ++j) {
            transposed.set(j, i) = get(i, j);
        }
    }

    return transposed;
}

/* RETURNS AN INVERSED MATRIX OBJECT */
template <size_t N>
Matrix<N, N> Matrix<N, N>::inverse() const {
    Matrix<N, N> A(*this);  // Copy of the matrix
    Matrix<N, N> I;         // Identity matrix
    I.identity();           // Initialize as identity

    for (size_t col = 0; col < N; col++) {
        // Step 1: Find pivot
        if (A.get(col, col) == 0) {
            // Swap with a row below
            bool swapped = false;
            for (size_t row = col + 1; row < N; row++) {
                if (A.get(row, col) != 0) {
                    A.swapRows(row, col);
                    I.swapRows(row, col);
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                throw std::runtime_error("Matrix is singular and cannot be inverted.");
            }
        }

        // Step 2: Normalize pivot row
        INT_DATATYPE pivot = A.get(col, col);
        for (size_t j = 0; j < N; j++) {
            A.set(col, j) /= pivot;
            I.set(col, j) /= pivot;
        }

        // Step 3: Zero out other column values
        for (size_t row = 0; row < N; row++) {
            if (row != col) {
                INT_DATATYPE factor = A.get(row, col);
                for (size_t j = 0; j < N; j++) {
                    A.set(row, j) -= factor * A.get(col, j);
                    I.set(row, j) -= factor * I.get(col, j);
                }
            }
        }
    }

    return I; // Inverse matrix
}
/* GAUSS-JORDAN HELPER FUNCTION */
template<size_t Rows, size_t Cols>
Matrix<N, N> Matrix<Rows, Cols>::gaussJordan() const {
    Matrix<N, N> A(*this);
    Matrix<N, N> I = Matrix<N, N>::identity();

    for (size_t col = 0; col < N; col++) {
        if(A(col, col) == 0) {
            bool swapped = false;
            for (size_t row = col + 1; row < N; row++) {
                if(A(row, col) != 0) {
                    std::swap(A[row], A[col]);
                    std::swap(I[row], I[col]);
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                throw std::runtime_error("Matrix is singular and cannot be inverted!")
            }
        }

        double pivot = A(col, col);
        for (size_t j = 0; j < N; j++) {
            A(col, j) /= pivot;
            A(col, j) /= pivot;
        }
        for (size_t row = 0; row < N; row++) {
            if (row != col) {
                double factor = A(row, col);
                for (size_t j = 0; j < N; j++) {
                    A(row, j) -= factor * A(col, j);
                    I(row, j) -= factor * I(col, j);
                }
            }
        }
    }
    return I;  
} 

/* RETURNS DETERMINANT OF MATRIX */ 
template <size_t Rows, size_t Cols>
INT_DATATYPE Matrix<Rows, Cols>::determinant() const {
    
    _SIZE_CHECK
    
    if constexpr (Rows == 1) {
        return get(0, 0);
    }
    else if constexpr (Rows == 2) {
        return get(0, 0) * get(1, 1) - get(0, 1) * get(1, 0);
    }
    else if constexpr (Rows == 3) {
        INT_DATATYPE A = get(0, 0) * (get(1, 1) * get(2, 2) - get(1, 2) * get(2, 1));
        INT_DATATYPE B = get(0, 1) * (get(1, 0) * get(2, 2) - get(1, 2) * get(2, 0));
        INT_DATATYPE C = get(0, 2) * (get(1, 0) * get(2, 1) - get(1, 1) * get(2, 0));
        return A - B + C;
    }
    else if constexpr (Rows <= 11) {
        INT_DATATYPE det = 0;
        for (size_t j = 0; j < Cols; ++j) {
            det += (j % 2 == 0 ? 1 : -1) * get(0, j) * minor_matrix(0, j).determinant();
        }
        return det;
    } else {
        INT_DATATYPE det = determinantLU(*this);
        return det;
    }
}

template <size_t Rows, size_t Cols>
Matrix<Rows - 1, Cols - 1> Matrix<Rows, Cols>::minor_matrix(size_t row, size_t col) const {
    static_assert(Rows > 1 && Cols > 1, "Minor matrix cannot be computed for matrices smaller than 2x2!");
    
    Matrix<Rows - 1, Cols - 1> minor;
    
    size_t minor_row = 0;
    for (size_t i = 0; i < Rows; ++i) {
        if (i == row) continue;
        
        size_t minor_col = 0;
        for (size_t j = 0; j < Cols; ++j) {
            if (j == col) continue;
            
            minor.set(minor_row, minor_col) = get(i, j);
            ++minor_col;
        }
        ++minor_row;
    }
    
    return minor;
}

template <size_t N>
INT_DATATYPE determinantLU(const Matrix<N, N>& mat) {
    Matrix<N, N> LU = mat;
    INT_DATATYPE det = 1;

    for (size_t k = 0; k < N; ++k) {
        size_t max_row = k;
        for (size_t i = k + 1; i < N; ++i) {
            if (abs(LU.get(i, k)) > abs(LU.get(max_row, k))) {
                max_row = i;
            }
        }
        if (max_row != k) {
            LU.swapRows(k, max_row);
            det *= -1;
        }

        if (LU.get(k, k) == 0) {
            return 0; 
        }

        det *= LU.get(k, k);

        for (size_t i = k + 1; i < N; ++i) {
            INT_DATATYPE factor = LU.get(i, k) / LU.get(k, k);
            for (size_t j = k + 1; j < N; ++j) {
                LU.set(i, j) -= factor * LU.get(k, j);
            }
        }
    }

    return det;
}

template <size_t Rows, size_t Cols>
void Matrix<Rows, Cols>::swapRows(size_t row1, size_t row2) {
    for (size_t j = 0; j < Cols; ++j) {
        std::swap(set(row1, j), set(row2, j));
    }
}


/* ADDITION OPERATOR */
template <size_t Rows, size_t Cols>
Matrix<Rows, Cols> Matrix<Rows, Cols>::operator+(const Matrix& other) const {
    Matrix result;

    for (size_t i = 0; i < Rows * Cols; ++i) {
        result.m_Data[i] = this->m_Data[i] + other.m_Data[i];
    }

    return result;
}

/* MATRIX SUBTRACTION */
template <size_t Rows, size_t Cols>
Matrix<Rows, Cols> Matrix<Rows, Cols>::operator-(const Matrix& other) const {

    Matrix result;

    for (size_t i = 0; i < Rows * Cols; i++) {
        result.m_Data[i] = this->m_Data[i] - other.m_Data[i];
    }

    return result;
}

/* MATRIX NEGATION */
template <size_t Rows, size_t Cols>
Matrix<Rows, Cols> Matrix<Rows, Cols>::operator-() const {
    Matrix result;

    for (size_t i = 0; i < Rows * Cols; i++) {
        result.m_Data[i] = this->m_Data[i] * -1;
    }

    return result;
}

/* MATRIX ON MATRIX MULTIPLICATION */
template <size_t Rows, size_t Cols>
template <size_t OtherCols>
Matrix<Rows, OtherCols> Matrix<Rows, Cols>::operator*(const Matrix<Cols, OtherCols>& other) const {
    static_assert(Cols == Rows, "Matrix dimensions are incompatible for multiplication");
    Matrix<Rows, OtherCols> result;

    for (size_t i = 0; i < Rows; ++i) {
        for (size_t j = 0; j < OtherCols; ++j) {
            result.set(i, j) = 0; // Initialize the result element to 0
            for (size_t k = 0; k < Cols; ++k) {
                result.set(i, j) += this->get(i, k) * other.get(k, j);
            }
        }
    }

    return result;
}


/* SCALAR MULTIPLICATION */
template <size_t Rows, size_t Cols>
Matrix<Rows, Cols> Matrix<Rows, Cols>::operator*(INT_DATATYPE scalar) const {
    Matrix result;

    for (size_t i = 0; i < Rows * Cols; ++i) {
        result.m_Data[i] = this->m_Data[i] * scalar;
    }

    return result;
}

template <size_t Rows, size_t Cols>
Matrix<Rows, Cols> operator*(INT_DATATYPE scalar, const Matrix<Rows, Cols>& matrix) {
    return matrix * scalar;
}


template <size_t Rows, size_t Cols>
size_t Matrix<Rows, Cols>::size() const {
    return Rows * Cols;
}

#endif // MATRIX_TPP
