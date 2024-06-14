/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/state/BoundSignaller.h>
#include <data/state/StateManager.h>

using namespace signaller;

BoundSignaller::BoundSignaller(unsigned int id, state::StateManager* const owner) : owner(owner), id(id) {}

void BoundSignaller::external_change() const {
    owner->externally_modified(id);
}

void BoundSignaller::internal_change() const {
    owner->internally_modified(id);
    owner->externally_modified(id);
}

unsigned int BoundSignaller::get_id() const {
    return id;
}