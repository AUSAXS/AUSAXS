#include <rigidbody/transform/SingleTransform.h>

using namespace rigidbody;

void SingleTransform::rotate(double rad, Constraint& constraint) {
    Body& body = constraint.get_body1();
    Vector3 r = constraint.get_atom1().coords - constraint.get_atom2().coords;
    Vector3<double> u1, u2, u3;
    std::tie(u1, u2, u3) = r.generate_basis();

    body.rotate(u1, rad);
}

void SingleTransform::translate(double length, Constraint& constraint) {
    Body& body = constraint.get_body1();
    Vector3 r = constraint.get_atom1().coords - constraint.get_atom2().coords;
    Vector3<double> u1, u2, u3;
    std::tie(u1, u2, u3) = r.generate_basis();

    body.translate(length*u1);
}