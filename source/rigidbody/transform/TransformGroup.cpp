#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <data/Body.h>

TransformGroup::TransformGroup(std::vector<Body*> bodies, std::vector<unsigned int> indices, std::shared_ptr<DistanceConstraint> target, Vector3<double> pivot) 
    : bodies(bodies), indices(indices), target(target), pivot(pivot) {}

TransformGroup::~TransformGroup() = default;
