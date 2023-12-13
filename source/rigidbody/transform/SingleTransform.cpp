#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/record/Atom.h>

using namespace rigidbody::transform;

SingleTransform::SingleTransform(RigidBody* rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(const Matrix<double>& M, const Vector3<double>& t, constraints::DistanceConstraint& constraint) {
    TransformGroup group({&constraint.get_body1()}, {constraint.ibody1}, constraint, constraint.get_atom2().coords);
    backup(group);
    rotate(M, group);
    translate(t, group);
}