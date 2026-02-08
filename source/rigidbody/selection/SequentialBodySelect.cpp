// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Random.h>

#include <random>

using namespace ausaxs::rigidbody::selection;

SequentialBodySelect::SequentialBodySelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialBodySelect::~SequentialBodySelect() = default;

std::pair<unsigned int, int> SequentialBodySelect::next() {
    unsigned int this_body = ibody;
    unsigned int N = rigidbody->constraints->discoverable_constraints.size();
    unsigned int M = rigidbody->constraints->get_body_constraints(this_body).size();
    ibody = (ibody + 1) % N;

    switch (M) {
        case 0: {
            return std::make_pair(this_body, -1);
        }
        case 1: {
            return std::make_pair(this_body, 0);
        }
        default: {
            std::uniform_int_distribution<int> distribution2(0, M-1);
            unsigned int iconstraint = distribution2(random::generator());

            return std::make_pair(this_body, iconstraint);
        }
    }
}
