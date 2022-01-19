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
        const vector<Constraint>& constraints = protein->constraints;

        std::unordered_map<size_t, size_t> group;
        std::for_each(protein->protein.bodies.begin(), protein->protein.bodies.end(), [&group] (const Body& body) {group[body.uid] = body.uid;});

        std::cout << "\n\nIDS: ";
        for (const auto& e : protein->protein.bodies) {
            std::cout << " " << e.uid;
        }
        std::cout << std::endl;

        for (size_t i = 0; i < constraints.size(); i++) {
            const Constraint& constraint = constraints[i];

            // if the constraint is the pivot
            if (constraint == pivot) {
                continue;
            }

            size_t id1 = group.at(constraint.body1->uid);
            size_t id2 = group.at(constraint.body2->uid);
            std::cout << "\tid1: " << id1 << ", id2: " << id2 << std::endl;
            if (group.at(id1) != group.at(id2)) {
                for (size_t i = 0; i < group.size(); i++) {
                    if (group.at(i) == id2) {
                        std::cout << "body " << i << " is part of group " << id2 << ", changing to group " << id1 << std::endl;
                        group[i] = id1;
                    }
                }

                // std::for_each(groups.begin(), groups.end(), [&groups, &id1, &id2] (const int& index) {if (groups[index] == id2) {groups[index] = id1;}});
            }

            std::for_each(protein->protein.bodies.begin(), protein->protein.bodies.end(), [&group] (const Body& body) {std::cout << group[body.uid];});
            std::cout << std::endl;
        }

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