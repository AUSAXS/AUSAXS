#include <rigidbody/selection/BodySelectStrategy.h>

#include <rigidbody/RigidBody.h>

using namespace rigidbody::selection;

BodySelectStrategy::BodySelectStrategy(const RigidBody* rigidbody) : rigidbody(rigidbody), N(rigidbody->body_size()) {}