#include <iostream>
#include <iomanip>
#include "matrix.hpp"

// Helper function to check if a matrix is approximately identity
template<typename T, size_t N>
bool isApproximatelyIdentity(const Matrix<T, N, N>& mat, T tolerance = 1e-6) {
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            T expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(mat.get(i, j) - expected) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

int main() {
    // Test case 1: Well-conditioned 3x3 matrix
    Matrix<double, 3, 3> mat1;
    mat1.set(0, 0, 4.0); mat1.set(0, 1, 7.0); mat1.set(0, 2, 2.0);
    mat1.set(1, 0, 3.0); mat1.set(1, 1, 5.0); mat1.set(1, 2, 6.0);
    mat1.set(2, 0, 8.0); mat1.set(2, 1, 1.0); mat1.set(2, 2, 9.0);

    // Test case 2: Poorly-conditioned matrix (Hilbert matrix 3x3)
    Matrix<double, 3, 3> mat2;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            mat2.set(i, j, 1.0 / (i + j + 1));
        }
    }

    std::cout << std::fixed << std::setprecision(15);
    
    std::cout << "Test Case 1 - Well-conditioned matrix:" << std::endl;
    std::cout << "Original matrix:" << std::endl;
    mat1.debug_print();
    
    auto inv1 = mat1.inverse();
    std::cout << "\nInverse:" << std::endl;
    inv1.debug_print();
    
    auto prod1 = mat1 * inv1;
    std::cout << "\nProduct (should be identity):" << std::endl;
    prod1.debug_print();
    
    std::cout << "\nIs approximately identity? " 
              << (isApproximatelyIdentity(prod1) ? "Yes" : "No") << std::endl;

    std::cout << "\nTest Case 2 - Hilbert matrix (poorly-conditioned):" << std::endl;
    std::cout << "Original matrix:" << std::endl;
    mat2.debug_print();
    
    auto inv2 = mat2.inverse();
    std::cout << "\nInverse:" << std::endl;
    inv2.debug_print();
    
    auto prod2 = mat2 * inv2;
    std::cout << "\nProduct (should be identity):" << std::endl;
    prod2.debug_print();
    
    std::cout << "\nIs approximately identity? " 
              << (isApproximatelyIdentity(prod2) ? "Yes" : "No") << std::endl;

    return 0;
}
