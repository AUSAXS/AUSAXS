#include <hist/detail/BodyTracker.h>
#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>

using namespace hist;

BodyTracker::BodyTracker(std::observer_ptr<const data::Molecule> protein) : body_size(protein->body_size()), statemanager(std::make_unique<state::StateManager>(body_size)) {}

BodyTracker::~BodyTracker() = default;

std::shared_ptr<signaller::Signaller> BodyTracker::get_probe(unsigned int i) {return statemanager->get_probe(i);}

void BodyTracker::signal_modified_hydration_layer() {
    statemanager->modified_hydration_layer();
}

std::observer_ptr<const state::StateManager> BodyTracker::get_state_manager() const {
    return std::make_observer(statemanager.get());
}

std::observer_ptr<state::StateManager> BodyTracker::get_state_manager() {
    return std::make_observer(statemanager.get());
}