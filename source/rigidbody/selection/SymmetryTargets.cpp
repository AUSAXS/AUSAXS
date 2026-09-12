// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SymmetryTargets.h>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/ReferenceSymmetry.h>

#include <algorithm>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::rigidbody::selection;

SymmetryTargets::SymmetryTargets(observer_ptr<const data::Molecule> molecule) : molecule(molecule) {
    assert(molecule != nullptr && "SymmetryTargets::SymmetryTargets: Molecule must not be null.");
}

void SymmetryTargets::refresh() const {
    if (!dirty) {return;}
    dirty = false;
    targets.clear();
    per_body.clear();

    for (int ibody = 0; ibody < molecule->size_body(); ++ibody) {
        const auto& body = molecule->get_body(ibody);
        for (int i = 0; i < body.size_symmetry(); ++i) {
            // is_optimizable sees through the ReferenceSymmetry wrapper to its base, and rejects the non-owning ReferenceSymmetryView
            if (!symmetry::is_optimizable(*body.symmetry().get(i))) {continue;}
            targets.emplace_back(Slot{.ibody=ibody, .isymmetry=i});
            per_body[ibody].emplace_back(i);
        }
    }
}

const std::vector<SymmetryTargets::Slot>& SymmetryTargets::all() const {
    refresh();
    return targets;
}

const std::vector<int>& SymmetryTargets::body_targets(int ibody) const {
    refresh();
    static const std::vector<int> none;
    auto it = per_body.find(ibody);
    return it == per_body.end() ? none : it->second;
}

std::optional<SymmetryTargets::Slot> SymmetryTargets::resolve(int ibody, int isymmetry) const {
    auto slot_exists = [&](int b, int i) {
        return b < molecule->size_body() && i < molecule->get_body(b).size_symmetry();
    };
    if (!slot_exists(ibody, isymmetry)) {return std::nullopt;}

    // a view carries no parameters of its own; the shared symmetry it stands for lives on the primary body, so redirect there. One hop suffices: a view always
    // refers to a ReferenceSymmetry, never to another view.
    const auto* sym = molecule->get_body(ibody).symmetry().get(isymmetry);
    if (const auto* view = dynamic_cast<const symmetry::ReferenceSymmetryView*>(sym)) {
        ibody = view->primary_body;
        isymmetry = view->symmetry_index;
        if (!slot_exists(ibody, isymmetry)) {return std::nullopt;}
    }

    const auto& drivable = body_targets(ibody);
    if (std::ranges::find(drivable, isymmetry) == drivable.end()) {return std::nullopt;}
    return Slot{.ibody=ibody, .isymmetry=isymmetry};
}
