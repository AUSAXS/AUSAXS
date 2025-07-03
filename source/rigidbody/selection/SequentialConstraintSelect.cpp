// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::selection;

SequentialConstraintSelect::SequentialConstraintSelect(observer_ptr<const RigidBody> rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialConstraintSelect::~SequentialConstraintSelect() = default;

std::pair<unsigned int, int> SequentialConstraintSelect::next() {
    unsigned int N = rigidbody->get_constraint_manager()->distance_constraints_map.size();
    unsigned int M = rigidbody->get_constraint_manager()->distance_constraints_map.at(ibody).size();

    if (iconstraint == M) {
        ibody = (ibody + 1) % N;
        iconstraint = 0;
    }

    return std::make_pair(ibody, iconstraint++);
}
