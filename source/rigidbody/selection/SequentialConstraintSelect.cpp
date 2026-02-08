// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::selection;

SequentialConstraintSelect::SequentialConstraintSelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialConstraintSelect::~SequentialConstraintSelect() = default;

std::pair<unsigned int, int> SequentialConstraintSelect::next() {
    unsigned int N = rigidbody->constraints->discoverable_constraints.size();
    unsigned int M = rigidbody->constraints->get_body_constraints(ibody).size();

    if (iconstraint == M) {
        ibody = (ibody + 1) % N;
        iconstraint = 0;
    }

    return std::make_pair(ibody, iconstraint++);
}
