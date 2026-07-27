// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/RandomBodySelect.h>
#include <utility/Random.h>

#include <random>

using namespace ausaxs::rigidbody::selection;

RandomBodySelect::RandomBodySelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

RandomBodySelect::~RandomBodySelect() = default;

BodySelectStrategy::Target RandomBodySelect::next(const ParameterMask& mask) {
    // under a symmetry-only mask the pose is frozen, so drawing a body would leave nothing to move whenever it holds only shared-symmetry views. Draw from
    // the drivable slots instead.
    if (symmetry_only(mask)) {return random_symmetry_target();}

    // rebuild the distribution each call against the live body count
    std::uniform_int_distribution<int> distribution(0, size_body()-1);
    unsigned int ibody = distribution(random::generator());
    return {ibody, random_constraint(ibody), -1};
}