#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Body.h>

TransformGroup::TransformGroup(std::vector<data::Body*> bodies, std::vector<unsigned int> indices, const rigidbody::DistanceConstraint& target, Vector3<double> pivot) 
    : bodies(bodies), indices(indices), target(target), pivot(pivot) {}

TransformGroup::~TransformGroup() = default;
