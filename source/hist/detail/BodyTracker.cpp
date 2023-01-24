#include <hist/detail/BodyTracker.h>
#include <data/Protein.h>

using namespace hist;

BodyTracker::BodyTracker(Protein* protein) : size(protein->bodies.size()), statemanager(size) {}

BodyTracker::~BodyTracker() = default;

std::shared_ptr<StateManager::BoundSignaller> BodyTracker::get_probe(unsigned int i) {return statemanager.get_probe(i);}

void BodyTracker::signal_modified_hydration_layer() {
    statemanager.modified_hydration_layer();
}

const StateManager& BodyTracker::get_state_manager() const {
    return statemanager;
}

StateManager& BodyTracker::get_state_manager() {
    return statemanager;
}