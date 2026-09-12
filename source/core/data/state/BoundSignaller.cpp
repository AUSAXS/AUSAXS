// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/state/BoundSignaller.h>

#include <data/state/StateManager.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::signaller;

BoundSignaller::BoundSignaller(int id, state::StateManager* const owner) : owner(owner), id(id) {}

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

void BoundSignaller::modified_hydration() const {
    owner->modified_hydration_layer();
}

void BoundSignaller::set_symmetry_size(int size) const {
    assert(static_cast<int>(owner->get_symmetry_modified_bodies()[id].size()) <= size && "BoundSignaller::set_symmetry_size: size is smaller than the current size");
    owner->get_symmetry_modified_bodies()[id].resize(size, true);
}

int BoundSignaller::get_symmetry_size() const {
    return static_cast<int>(owner->get_symmetry_modified_bodies()[id].size());
}

int BoundSignaller::get_id() const {
    return id;
}