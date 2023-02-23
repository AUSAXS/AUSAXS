#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody;

RigidTransform::RigidTransform(RigidBody* rigidbody) : TransformStrategy(rigidbody) {}

RigidTransform::~RigidTransform() = default;

void RigidTransform::rotate(const Matrix<double>& M, std::shared_ptr<Constraint> constraint) {
    auto group = get_connected(constraint);

    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (Body* body) {body->translate(-group.pivot);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&M] (Body* body) {body->rotate(M);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (Body* body) {body->translate(group.pivot);});
}

void RigidTransform::translate(const Vector3<double>& t, std::shared_ptr<Constraint> constraint) {
    auto group = get_connected(constraint);
    std::for_each(group.bodies.begin(), group.bodies.end(), [&t] (Body* body) {body->translate(t);});
}
