/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/state/BoundSignaller.h>
#include <data/state/StateManager.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::signaller;

BoundSignaller::BoundSignaller(unsigned int id, state::StateManager* const owner) : owner(owner), id(id) {}

void BoundSignaller::modified_external() const {
    owner->externally_modified(id);
}

void BoundSignaller::modified_internal() const {
    owner->internally_modified(id);
    owner->externally_modified(id);
}

void BoundSignaller::modified_symmetry(int i) const {
    owner->modified_symmetry(id, i);
}

void BoundSignaller::set_symmetry_size(std::size_t size) const {
    assert(owner->get_symmetry_modified_bodies()[id].size() <= size && "BoundSignaller::set_symmetry_size: size is smaller than the current size");
    owner->get_symmetry_modified_bodies()[id].resize(size, true);
}

unsigned int BoundSignaller::get_id() const {
    return id;
}