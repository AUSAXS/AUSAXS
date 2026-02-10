// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::state;

StateManager::StateManager(std::size_t size) 
    : _size(size), _externally_modified(size, true), _internally_modified(size, true), 
      _symmetry_modified(size, std::vector<bool>()), _modified_hydration(true), _modified(true)
{
    for (int i = 0; i < static_cast<int>(size); ++i) {
        probes.emplace_back(std::make_shared<signaller::BoundSignaller>(i, this));
    }
}

StateManager::StateManager(std::size_t size, const std::vector<std::size_t>& symmetry_sizes) 
    : _size(size), _externally_modified(size, true), _internally_modified(size, true), 
      _symmetry_modified(size), _modified_hydration(true), _modified(true)
{
    assert(symmetry_sizes.size() == size && "StateManager::StateManager: symmetry_sizes size mismatch");
    for (int i = 0; i < static_cast<int>(size); ++i) {
        _symmetry_modified[i] = std::vector<bool>(symmetry_sizes[i], false);
        probes.emplace_back(std::make_shared<signaller::BoundSignaller>(i, this));
    }
}

void StateManager::externally_modified_all() {
    _externally_modified = std::vector<bool>(size(), true);
    _modified = true;
}

void StateManager::internally_modified_all() {
    _internally_modified = std::vector<bool>(size(), true);
    _modified = true;
}

void StateManager::externally_modified(int i) {
    assert(i < static_cast<int>(size()) && "StateManager::externally_modified: index out of range");
    _externally_modified[i] = true;
    _modified = true;
}

void StateManager::internally_modified(int i) {
    assert(i < static_cast<int>(size()) && "StateManager::internally_modified: index out of range");
    _internally_modified[i] = true;
    _modified = true;
}

void StateManager::modified_hydration_layer() {
    _modified_hydration = true;
    _modified = true;
}

void StateManager::modified_symmetry(int i, int j) {
    assert(i < static_cast<int>(size()) && "StateManager::modified_symmetry: index out of range");
    assert(j < static_cast<int>(_symmetry_modified[i].size()) && "StateManager::modified_symmetry: index out of range");
    _symmetry_modified[i][j] = true;
    _modified = true;
}

void StateManager::reset_to_false() {
    _internally_modified = std::vector<bool>(size(), false);
    _externally_modified = std::vector<bool>(size(), false);
    _modified_hydration = false;
    _modified = false;
    for (auto& v : _symmetry_modified) {v = std::vector<bool>(v.size(), false);}
}

void StateManager::set_probe(int i, std::shared_ptr<signaller::Signaller> probe) {
    assert(i < static_cast<int>(probes.size()) && "StateManager::set_probe: index out of range");
    probes[i] = std::move(probe);
}

std::shared_ptr<signaller::Signaller> StateManager::get_probe(int i) {
    assert(i < static_cast<int>(probes.size()) && "StateManager::get_probe: index out of range");
    return probes[i];
}

std::vector<std::shared_ptr<signaller::Signaller>> StateManager::get_probes() {return probes;}

const std::vector<bool>& StateManager::get_externally_modified_bodies() const {return _externally_modified;}
std::vector<bool>& StateManager::get_externally_modified_bodies() {return _externally_modified;}

const std::vector<bool>& StateManager::get_internally_modified_bodies() const {return _internally_modified;}
std::vector<bool>& StateManager::get_internally_modified_bodies() {return _internally_modified;}

const std::vector<std::vector<bool>>& StateManager::get_symmetry_modified_bodies() const {return _symmetry_modified;}
std::vector<std::vector<bool>>& StateManager::get_symmetry_modified_bodies() {return _symmetry_modified;}

bool StateManager::is_modified() const {
    return _modified;
}

bool StateManager::is_externally_modified(int i) const {return _externally_modified[i];}

bool StateManager::is_internally_modified(int i) const {return _internally_modified[i];}

bool StateManager::is_modified_symmetry(int i, int j) const {return _symmetry_modified[i][j];}

bool StateManager::is_modified_hydration() const {return _modified_hydration;}

std::size_t StateManager::size() const {return _size;}