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
    vector<const Body*> get_connected(const Constraint& pivot) const {
        const vector<Constraint>& constraints = protein->constraints;       // easy access to the set of all constraints
        std::unordered_map<size_t, size_t> group;                           // a map of body uids to group ids
        std::unordered_map<size_t, std::list<const Body*>> group_members;   // a map of group ids to a list of its members

        // each body starts in its own group
        std::for_each(protein->protein.bodies.begin(), protein->protein.bodies.end(), [&group, &group_members] (const Body& body) {
            group[body.uid] = body.uid;
            group_members[body.uid] = {&body};
        });

        // iterate through all constraints
        for (size_t i = 0; i < constraints.size(); i++) {
            const Constraint& constraint = constraints[i];

            // if the constraint is the pivot, we skip it
            if (constraint == pivot) {
                continue;
            }

            // get the current group id of each body from the constraint
            size_t id1 = group.at(constraint.body1->uid);
            size_t id2 = group.at(constraint.body2->uid);

            // if they are not already in the same group, we merge their groups
            if (group.at(id1) != group.at(id2)) {
                // references to the two group lists
                std::list<const Body*>& members1 = group_members.at(id1);
                std::list<const Body*>& members2 = group_members.at(id2);

                // change all members of group2 to group1
                for (const auto& body : members2) {
                    group[body->uid] = id1;
                    members1.push_back(body);
                }
                group_members.erase(id2); // erase the old member list
            }
        }

        // get the id of the pivot, and return a vector of all bodies from the same group
        size_t id = group.at(pivot.body1->uid);
        vector<const Body*> connected;
        for(const auto& body : protein->protein.bodies) {
            if (group.at(body.uid) == id) {
                connected.push_back(&body);
            }
        }
        return connected;
    }

    /**
     * @brief Translate a body. 
     */
    void translate(const Vector3& v, Body& body) override {}
};