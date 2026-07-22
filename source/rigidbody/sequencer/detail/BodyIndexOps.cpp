// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>

#include <algorithm>

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

    std::sort(indices.begin(), indices.end());
    erase_indices(molecule.get_bodies(), indices);
    erase_indices(conformation.initial_conformation, indices);
    erase_indices(conformation.absolute_parameters.parameters, indices);

    owner->setup()._get_body_names().remove(indices);
}
