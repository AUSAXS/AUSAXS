// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/SymmetryTargets.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::rigidbody::selection;

SymmetryTargets::SymmetryTargets(observer_ptr<const data::Molecule> molecule) : molecule(molecule) {
    assert(molecule != nullptr && "SymmetryTargets::SymmetryTargets: Molecule must not be null.");
}

void SymmetryTargets::refresh() {
    if (!dirty) {return;}
    dirty = false;
    targets.clear();
    per_body.clear();

    for (unsigned int ibody = 0; ibody < molecule->size_body(); ++ibody) {
        const auto& body = molecule->get_body(ibody);
        for (unsigned int i = 0; i < body.size_symmetry(); ++i) {
            // is_optimizable sees through the ReferenceSymmetry wrapper to its base, and rejects the non-owning ReferenceSymmetryView
            if (!symmetry::is_optimizable(*body.symmetry().get(i))) {continue;}
            targets.emplace_back(Target{ibody, i});
            per_body[ibody].emplace_back(i);
        }
    }
}

const std::vector<SymmetryTargets::Target>& SymmetryTargets::all() {
    refresh();
    return targets;
}

const std::vector<unsigned int>& SymmetryTargets::body_targets(unsigned int ibody) {
    refresh();
    static const std::vector<unsigned int> none;
    auto it = per_body.find(ibody);
    return it == per_body.end() ? none : it->second;
}
