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
    void rotate(const Vector3& axis, const double rad, Body& body) override {
        std::for_each(protein->constraints.begin(), protein->constraints.end(), [] (const Constraint& constraint) {});

        body.rotate(axis, rad);
    }

    /**
     * @brief Get all bodies connected by constraints to the first body of the pivot. 
     *        If we have the four bodies A - B - C - D and pivot around the BC connection, this would return the group {AB}.
     *        TODO: redo this method in linear time, I know it must be possible. 
     */
    vector<Body> get_connected(const Constraint& pivot) const {
        const vector<Constraint>& constraints = protein->constraints;

        vector<size_t> groups(Body::uid_counter, -1);
        for (size_t i = 0; i < protein->protein.bodies.size(); i++) {
            groups[i] = i;
        }

        for (size_t i = 0; i < constraints.size(); i++) {
            const Constraint& constraint = constraints[i];

            // if the constraint is the pivot
            if (constraint == pivot) {
                continue;
            }

            size_t id1 = constraint.body1->uid;
            size_t id2 = constraint.body2->uid;
            if (groups[id1] != groups[id2]) {
                std::for_each(groups.begin(), groups.end(), [&groups, &id1, &id2] (const int& index) {if (groups[index] == id2) {groups[index] = id1;}});
            }
        }

        size_t id = groups[pivot.body1->uid];
        vector<Body> group;
        std::copy_if(protein->protein.bodies.begin(), protein->protein.bodies.end(), group.begin(), [&groups, &id] (const Body& body) {return groups[body.uid] == id;});
        return group;
    }

    /**
     * @brief Translate a body. 
     */
    void translate(const Vector3& v, Body& body) override {}
};