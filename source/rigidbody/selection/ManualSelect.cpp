// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/selection/SymmetryTargets.h>
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

BodySelectStrategy::Target ManualSelect::next(const ParameterMask& mask) {
    if (0 <= isymmetry) {return {ibody, -1, isymmetry};}

    // selecting the body itself leaves target_symmetry unset, so all of its symmetries move together - which is what "optimize this body" should mean, and
    // harmless for the undrivable views among them as long as at least one slot can move. If none can, and the mask has frozen the pose as well, the whole
    // run would perturb nothing at all: refuse rather than spend the optimizer's budget on it.
    if (symmetry_only(mask) && symmetry_candidates(ibody).empty()) {
        throw except::invalid_argument(
            "ManualSelect::next: The selected parameter mask optimizes symmetries only, but body " + std::to_string(ibody) + " declares none that can be "
            "driven. Note that a symmetry shared between several bodies is only drivable through the body owning it."
        );
    }
    return {ibody, random_constraint(ibody), -1};
}
