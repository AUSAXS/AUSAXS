#include <rigidbody/transform/SingleTransform.h>

using namespace rigidbody;

void SingleTransform::apply(const Matrix<double>& M, const Vector3<double>& t, std::shared_ptr<DistanceConstraint> constraint) {
    TransformGroup group({&constraint->get_body1()}, {constraint->ibody1}, constraint, constraint->get_atom1().coords);
    backup(group);
    rotate(M, group);
    translate(t, group);
}