// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>

#include <random>

using namespace ausaxs::rigidbody::selection;

SequentialBodySelect::SequentialBodySelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

SequentialBodySelect::~SequentialBodySelect() = default;

std::pair<unsigned int, int> SequentialBodySelect::next() {
    unsigned int this_body = ibody;
    unsigned int N = rigidbody->constraints->distance_constraints_map.size();
    unsigned int M = rigidbody->constraints->distance_constraints_map.at(this_body).size();
    ibody = (ibody + 1) % N;

    switch (M) {
        case 0: {
            return std::make_pair(this_body, -1);
        }
        case 1: {
            return std::make_pair(this_body, 0);
        }
        default: {
            std::random_device random;
            std::mt19937 generator2(random());
            std::uniform_int_distribution<int> distribution2(0, M-1);
            unsigned int iconstraint = distribution2(generator2);

            return std::make_pair(this_body, iconstraint);
        }
    }
}
