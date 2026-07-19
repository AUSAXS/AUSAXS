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

namespace {
    template<typename T>
    void erase_indices(std::vector<T>& vec, const std::vector<int>& sorted_indices) {
        if (sorted_indices.empty()) {return;}
        std::vector<T> kept;
        kept.reserve(vec.size() - sorted_indices.size());
        std::size_t next_removed = 0;
        for (std::size_t i = 0; i < vec.size(); ++i) {
            if (next_removed < sorted_indices.size() && sorted_indices[next_removed] == static_cast<int>(i)) {
                ++next_removed;
                continue;
            }
            kept.emplace_back(std::move(vec[i]));
        }
        vec = std::move(kept);
    }
}

void ausaxs::rigidbody::sequencer::detail::erase_bodies(observer_ptr<Sequencer> owner, std::vector<int> indices) {
    auto& molecule = *owner->_get_molecule();
    auto& conformation = *owner->_get_rigidbody()->conformation;
    auto& body_names = owner->setup()._get_body_names();

    std::sort(indices.begin(), indices.end());
    erase_indices(molecule.get_bodies(), indices);
    erase_indices(conformation.initial_conformation, indices);
    erase_indices(conformation.absolute_parameters.parameters, indices);

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
