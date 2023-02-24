#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody;

RigidTransform::RigidTransform(RigidBody* rigidbody) : TransformStrategy(rigidbody) {}

RigidTransform::~RigidTransform() = default;

void RigidTransform::apply(const Matrix<double>& M, const Vector3<double>& t, std::shared_ptr<Constraint> constraint) {
    auto group = get_connected(constraint);
    backup(group);

    // remove the bodies from the grid
    auto grid = rigidbody->get_grid();
    for (auto& body : group.bodies) {
        grid->remove(body);
    }

    rotate(M, group);
    translate(t, group);

    // add them back to the grid
    for (auto& body : group.bodies) {
        grid->add(body);
    }
}