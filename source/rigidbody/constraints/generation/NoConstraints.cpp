#include <rigidbody/constraints/generation/NoConstraints.h>
#include <rigidbody/constraints/DistanceConstraint.h>

using namespace rigidbody::constraints;

std::vector<DistanceConstraint> NoConstraints::generate() const {return {};}