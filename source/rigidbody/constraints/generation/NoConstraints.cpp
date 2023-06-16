#include <rigidbody/constraints/generation/NoConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>

using namespace rigidbody;

std::vector<rigidbody::DistanceConstraint> NoConstraints::generate() const {return {};}