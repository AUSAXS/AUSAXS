#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Body.h>

using namespace rigidbody::transform;

TransformGroup::TransformGroup(std::vector<data::Body*> bodies, std::vector<unsigned int> indices, const constraints::DistanceConstraint& target, Vector3<double> pivot) 
    : bodies(bodies), indices(indices), target(target), pivot(pivot) {}

TransformGroup::~TransformGroup() = default;