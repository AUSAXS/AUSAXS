#include <rigidbody/transform/SingleTransform.h>

using namespace rigidbody;

void SingleTransform::rotate(const Matrix<double>& M, std::shared_ptr<Constraint> constraint) {
    Body& body = constraint->get_body1();

    body.translate(-constraint->get_atom1().coords);
    body.rotate(M);
    body.translate(constraint->get_atom1().coords);
}

void SingleTransform::translate(const Vector3<double>& t, std::shared_ptr<Constraint> constraint) {
    Body& body = constraint->get_body1();
    body.translate(t);
}