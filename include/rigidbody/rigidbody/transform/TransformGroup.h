#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <math/Vector3.h>
#include <data/DataFwd.h>

#include <vector>

namespace ausaxs::rigidbody::transform {
    struct TransformGroup {
        TransformGroup(std::vector<data::Body*> bodies, std::vector<unsigned int> indices, const constraints::DistanceConstraint& target, Vector3<double> pivot);
        ~TransformGroup();
        std::vector<data::Body*> bodies;                // The bodies to transform.
        std::vector<unsigned int> indices;              // The indices of the bodies in the rigidbody.
        const constraints::DistanceConstraint& target;  // The constraint to transform along.
        Vector3<double> pivot;                          // The pivot point of the transformation.
    };
}