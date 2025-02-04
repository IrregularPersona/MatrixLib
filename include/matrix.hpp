#ifndef MATRIX_H
#define MATRIX_H

#include <array>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#define INT_DATATYPE int64_t
#define FLOAT_DATATYPE float
#define DOUBLE_DATATYPE double
// typedef __int128_t INT_DATATYPE;

template <typename T, size_t Rows, size_t Cols>
class Matrix {
private:
    std::array<T, Rows * Cols> m_Data;

public:
    Matrix() { std::fill(m_Data.begin(), m_Data.end(), 0); }

    void set(size_t i, size_t j, T value) { m_Data[i * Cols + j] = value; }
    T get(size_t i, size_t j) const { return m_Data[i * Cols + j]; }

    void identity() {
        static_assert(Rows == Cols, "Identity matrix must be square");
        std::fill(m_Data.begin(), m_Data.end(), 0);
        for (size_t i = 0; i < Rows; ++i) set(i, i, 1);
    }

    void debug_print() const {
        std::cout.precision(4);
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                std::cout << get(i, j) << " ";
            }
            std::cout << '\n';
        }
    }

    Matrix& init(T value) {
        std::fill(m_Data.begin(), m_Data.end(), value);
        return *this;
    }

    Matrix<T, Cols, Rows> transpose() const {
        Matrix<T, Cols, Rows> result;
        for (size_t i = 0; i < Rows; ++i)
            for (size_t j = 0; j < Cols; ++j)
                result.set(j, i, get(i, j));
        return result;
    }

    Matrix<T, Rows, Cols> inverse() const {
        static_assert(Rows == Cols, "Inverse requires square matrix");
        const size_t N = Rows;

        Matrix<T, N, 2 * N> augmented;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                augmented.set(i, j, get(i, j)); // Copy original matrix
                augmented.set(i, j + N, (i == j) ? 1 : 0); // Identity matrix
            }
        }

        for (size_t pivot = 0; pivot < N; ++pivot) {
            size_t max_row = pivot;
            for (size_t i = pivot + 1; i < N; ++i) {
                if (std::abs(augmented.get(i, pivot)) > std::abs(augmented.get(max_row, pivot))) {
                    max_row = i;
                }
            }

            if (max_row != pivot) {
                augmented.swapRows(pivot, max_row);
            }

            if (augmented.get(pivot, pivot) == 0) {
                throw std::runtime_error("Matrix is singular and cannot be inverted");
            }

            T pivot_value = augmented.get(pivot, pivot);
            for (size_t j = pivot; j < 2 * N; ++j) {
                augmented.set(pivot, j, augmented.get(pivot, j) / pivot_value);
            }

            for (size_t i = 0; i < N; ++i) {
                if (i != pivot) {
                    T factor = augmented.get(i, pivot);
                    for (size_t j = pivot; j < 2 * N; ++j) {
                        augmented.set(i, j, augmented.get(i, j) - factor * augmented.get(pivot, j));
                    }
                }
            }
        }

        // Extract the inverse matrix from the augmented matrix
        Matrix<T, N, N> result;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                result.set(i, j, augmented.get(i, j + N));
            }
        }

        return result;
    }

    T determinant() const {
        static_assert(Rows == Cols, "Determinant requires square matrix");
        Matrix<T, Rows, Cols> LU = *this;
        T det = 1;

        for (size_t k = 0; k < Rows; ++k) {
            size_t max_row = k;
            for (size_t i = k + 1; i < Rows; ++i) {
                if (std::abs(LU.get(i, k)) > std::abs(LU.get(max_row, k))) {
                    max_row = i;
                }
            }
            
            if (max_row != k) {
                LU.swapRows(k, max_row);
                det *= -1;
            }

            if (LU.get(k, k) == 0) return 0;

            det *= LU.get(k, k);
            
            for (size_t i = k + 1; i < Rows; ++i) {
                T factor = LU.get(i, k) / LU.get(k, k);
                for (size_t j = k + 1; j < Cols; ++j) {
                    LU.set(i, j, LU.get(i, j) - factor * LU.get(k, j));
                }
            }
        }

        return det;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix result;
        for (size_t i = 0; i < Rows * Cols; ++i)
            result.m_Data[i] = m_Data[i] + other.m_Data[i];
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result;
        for (size_t i = 0; i < Rows * Cols; ++i)
            result.m_Data[i] = m_Data[i] - other.m_Data[i];
        return result;
    }

    Matrix operator-() const {
        Matrix result;
        for (size_t i = 0; i < Rows * Cols; ++i)
            result.m_Data[i] = -m_Data[i];
        return result;
    }
    
    Matrix operator*(T scalar) const {
        Matrix result;
        for (size_t i = 0; i < Rows * Cols; ++i)
            result.m_Data[i] = m_Data[i] * scalar;
        return result;
    }

    template<size_t OtherCols>
    Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols>& other) const {
        Matrix<T, Rows, OtherCols> result;
        result.init(0);
    
        for (size_t i = 0; i < Rows; i++) {
            for (size_t j = 0; j < OtherCols; j++) {
                T sum = 0;
                T c = 0;  // compensation term
                for (size_t k = 0; k < Cols; k++) {
                    T product = get(i, k) * other.get(k, j);
                    T y = product - c;
                    T t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                result.set(i, j, sum);
            }
        }
    
        return result;
    }
    
    void swapRows(size_t row1, size_t row2) {
        for (size_t j = 0; j < Cols; ++j) {
            std::swap(m_Data[row1 * Cols + j], m_Data[row2 * Cols + j]);
        }
    }

    size_t size() const { return Rows * Cols; }
};

template <typename T, size_t Rows, size_t Cols>
Matrix<T, Rows, Cols> operator*(T scalar, const Matrix<T, Rows, Cols>& matrix) {
    return matrix * scalar;
}

#endif // MATRIX_H
