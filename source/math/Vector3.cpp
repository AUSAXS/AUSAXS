#include "math/Vector3.h"
#include "math/Matrix.h"

void Vector3::rotate(const Matrix& matrix) {
    *this = matrix*(*this);
}

void Vector3::rotate(const Vector3& axis, double angle) {
    Matrix R = Matrix::rotation_matrix(axis, angle);
    rotate(R);
}
