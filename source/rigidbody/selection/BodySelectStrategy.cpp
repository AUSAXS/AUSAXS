// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/selection/ParameterMaskStrategy.h>

#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::selection;

BodySelectStrategy::BodySelectStrategy(observer_ptr<const Rigidbody> rigidbody)
    : rigidbody(rigidbody), N(rigidbody->molecule.size_body()),
      mask_strategy(std::make_unique<AllMaskStrategy>())
{}

BodySelectStrategy::SelectionResult BodySelectStrategy::next_mask() {
    auto [ibody, iconstraint] = next();
    return {ibody, iconstraint, mask_strategy->next()};
}

void BodySelectStrategy::set_mask_strategy(std::unique_ptr<ParameterMaskStrategy> strategy) {
    mask_strategy = std::move(strategy);
}