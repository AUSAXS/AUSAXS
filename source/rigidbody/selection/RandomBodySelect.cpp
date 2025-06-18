// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Exceptions.h>
#include <utility/Random.h>

using namespace ausaxs::rigidbody::selection;

RandomBodySelect::RandomBodySelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {
    distribution = std::uniform_int_distribution<int>(0, N-1);
}

RandomBodySelect::~RandomBodySelect() = default;

std::pair<unsigned int, int> RandomBodySelect::next() {
    unsigned int ibody = distribution(random::generator());

    unsigned int N = rigidbody->constraints->distance_constraints_map.at(ibody).size();
    switch (N) {
        case 0: {
            return std::make_pair(ibody, -1);
        }
        case 1: {
            return std::make_pair(ibody, 0);
        }
        default: {
            std::uniform_int_distribution<int> distribution2(0, N-1);
            unsigned int iconstraint = distribution2(random::generator());

            return std::make_pair(ibody, iconstraint);
        }
    }
}