#pragma once

#include <math/Vector3.h>
#include <data/DataFwd.h>

#include <vector>

namespace rigidbody {class DistanceConstraint;}
struct TransformGroup {
    TransformGroup(std::vector<data::Body*> bodies, std::vector<unsigned int> indices, const rigidbody::DistanceConstraint& target, Vector3<double> pivot);
    ~TransformGroup();
    std::vector<data::Body*> bodies;                // The bodies to transform.
    std::vector<unsigned int> indices;              // The indices of the bodies in the rigidbody.
    const rigidbody::DistanceConstraint& target;    // The constraint to transform along.
    Vector3<double> pivot;                          // The pivot point of the transformation.
};