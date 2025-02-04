#include "matrix.hpp"

int main() {
 Matrix<float, 3, 3> mat1;
 mat1.init(3);
 mat1.debug_print();
 return 0;
}
