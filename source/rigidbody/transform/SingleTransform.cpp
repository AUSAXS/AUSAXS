#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Atom.h>

using namespace rigidbody;

SingleTransform::SingleTransform(RigidBody* rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(const Matrix<double>& M, const Vector3<double>& t, std::shared_ptr<DistanceConstraint> constraint) {
    TransformGroup group({&constraint->get_body1()}, {constraint->ibody1}, constraint, constraint->get_atom1().coords);
    backup(group);
    rotate(M, group);
    translate(t, group);
}