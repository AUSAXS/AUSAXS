// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Random.h>

#include <random>

using namespace ausaxs::rigidbody::selection;

ManualSelect::ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody, bool use_constraints)
    : BodySelectStrategy(rigidbody), ibody(ibody), use_constraints(use_constraints)
{}

ManualSelect::~ManualSelect() = default;

std::pair<unsigned int, int> ManualSelect::next() {
    if (!use_constraints) {return std::make_pair(ibody, -1);}

    unsigned int N = rigidbody->constraints->get_body_constraints(ibody).size();
    switch (N) {
        case 0: {
            return std::make_pair(ibody, -1);
        }
        case 1: {
            return std::make_pair(ibody, 0);
        }
        default: {
            std::uniform_int_distribution<int> distribution(0, N-1);
            unsigned int iconstraint = distribution(random::generator());
            return std::make_pair(ibody, iconstraint);
        }
    }
}
