#include <hist/detail/BodyTracker.h>
#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>

using namespace hist;

BodyTracker::BodyTracker(Protein* protein) : size(protein->body_size()), statemanager(std::make_unique<StateManager>(size)) {}

BodyTracker::~BodyTracker() = default;

std::shared_ptr<signaller::Signaller> BodyTracker::get_probe(unsigned int i) {return statemanager->get_probe(i);}

void BodyTracker::signal_modified_hydration_layer() {
    statemanager->modified_hydration_layer();
}

const StateManager& BodyTracker::get_state_manager() const {
    return *statemanager;
}

StateManager& BodyTracker::get_state_manager() {
    return *statemanager;
}