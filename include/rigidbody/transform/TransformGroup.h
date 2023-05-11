#pragma once

#include <math/Vector3.h>

#include <vector>

class Body;
namespace rigidbody {class DistanceConstraint;}
struct TransformGroup {
    TransformGroup(std::vector<Body*> bodies, std::vector<unsigned int> indices, std::shared_ptr<rigidbody::DistanceConstraint> target, Vector3<double> pivot);
    ~TransformGroup();
    std::vector<Body*> bodies;                              // The bodies to transform.
    std::vector<unsigned int> indices;                      // The indices of the bodies in the rigidbody.
    std::shared_ptr<rigidbody::DistanceConstraint> target;  // The constraint to transform along.
    Vector3<double> pivot;                                  // The pivot point of the transformation.
};