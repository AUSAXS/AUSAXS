#pragma once

#include <array>
#include <algorithm>

#include "data/Protein.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/TransformationStrategy.h"

/**
 * @brief \class RigidTransform. 
 * 
 * With this transformation strategy, everything connected to the target of the transformation will be transformed as well. 
 */
class RigidTransform : public TransformationStrategy {
  public:
    /**
     * @brief Construtor. 
     */
    RigidTransform(const RigidBody* protein) : TransformationStrategy(protein) {}

    /**
     * @brief Destructor.
     */
    ~RigidTransform() override = default;

    /**
     * @brief Rotate a body. 
     */
    void rotate(const Vector3& axis, const double rad, Constraint& constraint) override {
        vector<Body*> bodies = get_connected(constraint);
    }

    /**
     * @brief Translate a body. 
     */
    void translate(const Vector3& v, Constraint& constraint) override {
        vector<Body*> bodies = get_connected(constraint);
        std::for_each(bodies.begin(), bodies.end(), [&v] (Body* const body) {body->translate(v);});
    }
};