// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Exceptions.h>
#include <utility/Random.h>

#include <random>

using namespace ausaxs::rigidbody::selection;

RandomBodySelect::RandomBodySelect(observer_ptr<const Rigidbody> rigidbody) : BodySelectStrategy(rigidbody) {}

RandomBodySelect::~RandomBodySelect() = default;

BodySelectStrategy::Target RandomBodySelect::next(const ParameterMask& mask) {
    // under a symmetry-only mask the pose is frozen, so drawing a body would leave nothing to move whenever it holds only shared-symmetry views. Draw from
    // the drivable slots instead; a symmetry move displaces none of the host body's own atoms, so there is no constraint group to drag along.
    if (symmetry_only(mask)) {
        const auto& candidates = symmetry_candidates();
        std::uniform_int_distribution<std::size_t> distribution(0, candidates.size()-1);
        const auto& target = candidates[distribution(random::generator())];
        return {target.ibody, -1, static_cast<int>(target.isymmetry)};
    }

    // rebuild the distribution each call against the live body count
    std::uniform_int_distribution<int> distribution(0, size_body()-1);
    unsigned int ibody = distribution(random::generator());
    return {ibody, random_constraint(ibody), -1};
}