// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/selection/ParameterMaskStrategy.h>
#include <rigidbody/selection/SymmetryTargets.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Exceptions.h>
#include <utility/Random.h>

#include <algorithm>
#include <cassert>
#include <random>

using namespace ausaxs::rigidbody::selection;

BodySelectStrategy::BodySelectStrategy(observer_ptr<const Rigidbody> rigidbody)
    : rigidbody(rigidbody),
      mask_strategy(std::make_unique<AllMaskStrategy>())
{}

unsigned int BodySelectStrategy::size_body() const {
    return rigidbody->molecule.size_body();
}

bool BodySelectStrategy::symmetry_only(const ParameterMask& mask) {
    return !mask.real_translation && !mask.real_rotation;
}

const std::vector<SymmetryTargets::Target>& BodySelectStrategy::symmetry_candidates() const {
    return rigidbody->symmetry_targets->all();
}

int BodySelectStrategy::random_constraint(unsigned int ibody) const {
    unsigned int N = rigidbody->constraints->get_body_constraints(ibody).size();
    if (N == 0) {return -1;}
    if (N == 1) {return 0;}
    std::uniform_int_distribution<int> distribution(0, N-1);
    return distribution(random::generator());
}

BodySelectStrategy::SelectionResult BodySelectStrategy::next_mask() {
    auto mask = mask_strategy->next();

    // a symmetry-only step can only accomplish something if there is a slot to drive; catching it here rather than letting every strategy discover an empty
    // pool on its own keeps the diagnostic in one place, and the alternative - quietly doing a pose move instead - would disobey the script
    if (symmetry_only(mask) && symmetry_candidates().empty()) {
        throw except::invalid_argument(
            "BodySelectStrategy::next_mask: The selected parameter mask optimizes symmetries only, but the molecule declares none that can be driven. "
            "Note that a symmetry shared between several bodies is only drivable through the body owning it."
        );
    }

    auto target = next(mask);

    // a target naming an undrivable slot would spend the whole iteration - hydration, histogram, fit - on a move that cannot change anything, and would do so
    // silently. Enforce here rather than in each strategy, so a new one cannot get it wrong unnoticed.
    assert(
        (target.isymmetry < 0 || [&] {
            const auto& drivable = rigidbody->symmetry_targets->body_targets(target.ibody);
            return std::find(drivable.begin(), drivable.end(), static_cast<unsigned int>(target.isymmetry)) != drivable.end();
        }())
        && "BodySelectStrategy::next_mask: the strategy selected a symmetry that cannot be driven."
    );

    // the selector is the authority on which symmetry moves; the mask merely carries that decision to ParameterMask::apply
    if (0 <= target.isymmetry) {mask.target_symmetry = static_cast<unsigned int>(target.isymmetry);}
    return {target.ibody, target.iconstraint, target.isymmetry, mask};
}

void BodySelectStrategy::set_mask_strategy(std::unique_ptr<ParameterMaskStrategy> strategy) {
    mask_strategy = std::move(strategy);
}
