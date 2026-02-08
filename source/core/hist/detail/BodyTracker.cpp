// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/BodyTracker.h>
#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::hist;

BodyTracker::BodyTracker(observer_ptr<const data::Molecule> protein) : body_size(protein->size_body()) {
    std::vector<std::size_t> symmetry_sizes(body_size);
    for (int i = 0; i < static_cast<int>(body_size); ++i) {
        symmetry_sizes[i] = protein->get_body(i).size_symmetry();
    }
    statemanager = std::make_unique<state::StateManager>(body_size, symmetry_sizes);
}

BodyTracker::~BodyTracker() = default;

std::shared_ptr<signaller::Signaller> BodyTracker::get_probe(unsigned int i) {return statemanager->get_probe(i);}

void BodyTracker::signal_modified_hydration_layer() {
    statemanager->modified_hydration_layer();
}

observer_ptr<const state::StateManager> BodyTracker::get_state_manager() const {
    return statemanager.get();
}

observer_ptr<state::StateManager> BodyTracker::get_state_manager() {
    return statemanager.get();
}