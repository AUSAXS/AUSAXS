// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/selection/SymmetryTargets.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <utility/Exceptions.h>

#include <algorithm>

using namespace ausaxs::rigidbody::selection;

ManualSelect::ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody) : BodySelectStrategy(rigidbody), ibody(ibody) {}

ManualSelect::ManualSelect(observer_ptr<const Rigidbody> rigidbody, unsigned int ibody, unsigned int isymmetry)
    : BodySelectStrategy(rigidbody), ibody(ibody), isymmetry(static_cast<int>(isymmetry))
{
    // reject an undrivable slot here rather than letting the optimizer spend its whole budget on a move that cannot change anything
    const auto& drivable = rigidbody->symmetry_targets->body_targets(ibody);
    if (std::find(drivable.begin(), drivable.end(), isymmetry) == drivable.end()) {
        throw except::invalid_argument(
            "ManualSelect::ManualSelect: Symmetry " + std::to_string(isymmetry) + " of body " + std::to_string(ibody) + " cannot be optimized. "
            "A symmetry shared between several bodies is only drivable through the body owning it."
        );
    }
}

ManualSelect::~ManualSelect() = default;

BodySelectStrategy::Target ManualSelect::next(const ParameterMask&) {
    if (0 <= isymmetry) {return {ibody, -1, isymmetry};}
    return {ibody, random_constraint(ibody), -1};
}
