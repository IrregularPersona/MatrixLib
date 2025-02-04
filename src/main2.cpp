#include "matrix.hpp"

int main() {
    Matrix<int, 2, 2> m1;
    Matrix<int, 2, 2> m2;

    // 1 5
    // 2 4

    // 10 7 
    // 12 3
    
    m1.set(0, 0, 1); m1.set(0, 1, 2); m1.set(1, 0, 5); m1.set(1, 1, 4);
    m2.set(0, 0, 10); m2.set(0, 1, 12); m2.set(1, 0, 7); m2.set(1, 1, 3);

    auto res = m1 * m2;

    res.debug_print();

    return 0;
    
}
