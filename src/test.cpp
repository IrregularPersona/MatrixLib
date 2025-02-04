#include <iostream>

#include "matrix.hpp"



int main() {

    Matrix<FLOAT_DATATYPE, 3, 3> mat;
    mat.set(0, 0, 4); mat.set(0, 1, 7); mat.set(0, 2, 2);
    mat.set(1, 0, 3); mat.set(1, 1, 5); mat.set(1, 2, 6);
    mat.set(2, 0, 8); mat.set(2, 1, 1); mat.set(2, 2, 9);

    std::cout << "Original Matrix:\n";

    mat.debug_print();
    try {
        auto inv = mat.inverse();
        std::cout << "\nInverse Matrix:\n";
        inv.debug_print();
        // Verify by multiplying: mat * inv should be the identity matrix

        auto product = mat * inv;

        std::cout << "\nProduct (should be identity):\n";

        product.debug_print();
        std::cout << "\nDetailed multiplication for element (1,0):\n";
        float sum = 0;
        sum += mat.get(1,0) * inv.get(0,0);
        std::cout << "3 * 0.1703 = " << sum << "\n";
        sum += mat.get(1,1) * inv.get(1,0);
        std::cout << "+ 5 * 0.0917 = " << sum << "\n";
        sum += mat.get(1,2) * inv.get(2,0);
        std::cout << "+ 6 * -0.1616 = " << sum << " (final)\n";
    } catch (const std::runtime_error& e) {

        std::cerr << "Error: " << e.what() << "\n";

    }


    return 0;

}

