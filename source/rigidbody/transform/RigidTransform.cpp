#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody;

RigidTransform::RigidTransform(RigidBody* rigidbody) : TransformStrategy(rigidbody) {}

RigidTransform::~RigidTransform() = default;

void RigidTransform::rotate(double rad, Constraint& constraint) {
    std::vector<Body*> bodies = get_connected(constraint);

    Vector3 r = constraint.get_atom1().coords - constraint.get_atom2().coords;
    Vector3<double> u1, u2, u3;
    std::tie(u1, u2, u3) = r.generate_basis();

    std::for_each(bodies.begin(), bodies.end(), [&u1, &rad] (Body* body) {body->rotate(u1, rad);});
}

void RigidTransform::translate(double length, Constraint& constraint) {
    std::vector<Body*> bodies = get_connected(constraint);

    Vector3 r = constraint.get_atom1().coords - constraint.get_atom2().coords;
    Vector3<double> u1, u2, u3;
    std::tie(u1, u2, u3) = r.generate_basis();

    std::for_each(bodies.begin(), bodies.end(), [&u1, &length] (Body* body) {body->translate(length*u1);});
}
