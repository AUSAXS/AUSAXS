#include <hist/detail/BodyTracker.h>
#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>

using namespace hist;

BodyTracker::BodyTracker(const data::Molecule* const protein) : body_size(protein->body_size()), statemanager(std::make_unique<state::StateManager>(body_size)) {}

BodyTracker::~BodyTracker() = default;

std::shared_ptr<signaller::Signaller> BodyTracker::get_probe(unsigned int i) {return statemanager->get_probe(i);}

void BodyTracker::signal_modified_hydration_layer() {
    statemanager->modified_hydration_layer();
}

const state::StateManager* BodyTracker::get_state_manager() const {
    return statemanager.get();
}

state::StateManager* BodyTracker::get_state_manager() {
    return statemanager.get();
}