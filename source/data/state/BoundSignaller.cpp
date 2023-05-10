#include <data/state/BoundSignaller.h>
#include <data/state/StateManager.h>

using namespace signaller;

BoundSignaller::BoundSignaller(unsigned int id, StateManager* const owner) : owner(owner), id(id) {}

void BoundSignaller::external_change() const {
    owner->externally_modified(id);
}

void BoundSignaller::internal_change() const {
    owner->internally_modified(id);
    owner->externally_modified(id);
}