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

const std::vector<SymmetryTargets::Slot>& BodySelectStrategy::symmetry_candidates() const {
    return rigidbody->symmetry_targets->all();
}

const std::vector<unsigned int>& BodySelectStrategy::symmetry_candidates(unsigned int ibody) const {
    return rigidbody->symmetry_targets->body_targets(ibody);
}

// the two helpers below wrap a slot as {ibody, -1, isymmetry}: driving a symmetry leaves the host body's own atoms where they are, so the move has no rigid
// group to propagate to and thus never carries a constraint. The empty-pool case is an assert rather than a throw because next_mask already rejects it with a
// user-facing diagnostic before any strategy is asked for a target; reaching it here means a caller went around that check.

BodySelectStrategy::Target BodySelectStrategy::random_symmetry_target() const {
    const auto& candidates = symmetry_candidates();
    assert(!candidates.empty() && "BodySelectStrategy::random_symmetry_target: no drivable symmetry to draw from.");
    std::uniform_int_distribution<std::size_t> distribution(0, candidates.size()-1);
    const auto& slot = candidates[distribution(random::generator())];
    return {slot.ibody, -1, static_cast<int>(slot.isymmetry)};
}

BodySelectStrategy::Target BodySelectStrategy::next_symmetry_target(std::size_t& cursor) const {
    const auto& candidates = symmetry_candidates();
    assert(!candidates.empty() && "BodySelectStrategy::next_symmetry_target: no drivable symmetry to step through.");

    // wrap on entry rather than after the increment: the pool can shrink between calls, so the cursor left behind last time may already be out of range
    cursor %= candidates.size();
    const auto& slot = candidates[cursor++];
    return {slot.ibody, -1, static_cast<int>(slot.isymmetry)};
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
            const auto& drivable = symmetry_candidates(target.ibody);
            return std::find(drivable.begin(), drivable.end(), static_cast<unsigned int>(target.isymmetry)) != drivable.end();
        }())
        && "BodySelectStrategy::next_mask: the strategy selected a symmetry that cannot be driven."
    );

    // the selector is the authority on which symmetry moves; the mask merely carries that decision to ParameterMask::apply
    if (0 <= target.isymmetry) {mask.target_symmetry = static_cast<unsigned int>(target.isymmetry);}
    return {target.ibody, target.iconstraint, std::move(mask)};
}

void BodySelectStrategy::set_mask_strategy(std::unique_ptr<ParameterMaskStrategy> strategy) {
    mask_strategy = std::move(strategy);
}
