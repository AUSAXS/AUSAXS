// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/selection/SymmetryTargets.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Exceptions.h>

using namespace ausaxs::rigidbody::selection;

ManualSelect::ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody) : BodySelectStrategy(rigidbody), ibody(ibody) {}

ManualSelect::ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody, unsigned int isymmetry) : BodySelectStrategy(rigidbody) {
    // a body participating in a shared symmetry declares it as a view onto the owner's copy. Naming that slot means the shared symmetry, so follow the link
    // rather than refusing the only name the user has for this body. Anything genuinely undrivable is rejected here, before the optimizer spends its whole
    // budget on a move that cannot change anything.
    auto target = rigidbody->symmetry_targets->resolve(ibody, isymmetry);
    if (!target.has_value()) {
        throw except::invalid_argument(
            "ManualSelect::ManualSelect: Symmetry " + std::to_string(isymmetry) + " of body " + std::to_string(ibody) + " cannot be optimized."
        );
    }
    this->ibody = target->ibody;
    this->isymmetry = static_cast<int>(target->isymmetry);
}

ManualSelect::~ManualSelect() = default;

BodySelectStrategy::Target ManualSelect::next(const ParameterMask&) {
    if (0 <= isymmetry) {return {ibody, -1, isymmetry};}
    return {ibody, random_constraint(ibody), -1};
}
