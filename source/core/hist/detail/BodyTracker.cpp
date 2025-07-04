// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/BodyTracker.h>
#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>
#include <data/Molecule.h>

using namespace ausaxs;
using namespace ausaxs::hist;

BodyTracker::BodyTracker(observer_ptr<const data::Molecule> protein) : body_size(protein->size_body()), statemanager(std::make_unique<state::StateManager>(body_size)) {}

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