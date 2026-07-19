// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>

#include <algorithm>
#include <iterator>

using namespace ausaxs::rigidbody::sequencer;

void ausaxs::rigidbody::sequencer::detail::erase_bodies(observer_ptr<Sequencer> owner, std::vector<int> indices) {
    auto& molecule = *owner->_get_molecule();
    auto& conformation = *owner->_get_rigidbody()->conformation;
    auto& body_names = owner->setup()._get_body_names();

    // erase highest index first so earlier indices stay valid while erasing
    std::sort(indices.begin(), indices.end());
    for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
        molecule.get_bodies().erase(molecule.get_bodies().begin() + *it);
        conformation.initial_conformation.erase(conformation.initial_conformation.begin() + *it);
        conformation.absolute_parameters.parameters.erase(conformation.absolute_parameters.parameters.begin() + *it);
    }

    // drop the erased bodies' names and reindex the survivors to their shifted body indices
    for (auto it = body_names.begin(); it != body_names.end(); ) {
        auto sel = detail::from_index(it->second);
        auto pos = std::lower_bound(indices.begin(), indices.end(), sel.body);
        if (pos != indices.end() && *pos == sel.body) {
            it = body_names.erase(it);
            continue;
        }
        auto shift = std::distance(indices.begin(), pos);
        if (shift != 0) {it->second = detail::to_index(sel.body - static_cast<int>(shift), sel.symmetry, sel.replica);}
        ++it;
    }
}
