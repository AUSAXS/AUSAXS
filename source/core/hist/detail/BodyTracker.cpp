/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/detail/BodyTracker.h>
#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>

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