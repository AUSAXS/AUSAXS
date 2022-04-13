#include <data/StateManager.h>

#include <memory>
#include <iostream>

StateManager::BoundSignaller::BoundSignaller(unsigned int id, StateManager* const owner) : owner(owner), id(id) {}

void StateManager::BoundSignaller::external_change() const {
    owner->externally_modified(id);
}

void StateManager::BoundSignaller::internal_change() const {
    owner->internally_modified(id);
    owner->externally_modified(id);
}

void StateManager::UnboundSignaller::external_change() const {}

void StateManager::UnboundSignaller::internal_change() const {}

StateManager::StateManager(unsigned int size) : size(size), _externally_modified(size, true), _internally_modified(size, true), _modified_hydration(true) {
    for (unsigned int i = 0; i < size; i++) {
        probes.push_back(std::make_shared<BoundSignaller>(i, this));
    }
}

void StateManager::externally_modified_all() {
    _externally_modified = std::vector<bool>(size, true);
}

void StateManager::internally_modified_all() {
    _internally_modified = std::vector<bool>(size, true);
}

void StateManager::externally_modified(const int i) {
    _externally_modified[i] = true;
}

void StateManager::internally_modified(const int i) {
    _internally_modified[i] = true;
}

void StateManager::modified_hydration_layer() {
    _modified_hydration = true;
}

void StateManager::reset() {
    _internally_modified = std::vector<bool>(size, false);
    _externally_modified = std::vector<bool>(size, false);
    _modified_hydration = false;
}

std::shared_ptr<StateManager::BoundSignaller> StateManager::get_probe(unsigned int i) {return probes[i];}

std::vector<bool> StateManager::get_externally_modified_bodies() const {return _externally_modified;}

std::vector<bool> StateManager::get_internally_modified_bodies() const {return _internally_modified;}

bool StateManager::is_externally_modified(unsigned int i) {return _externally_modified[i];}

bool StateManager::is_internally_modified(unsigned int i) {return _internally_modified[i];}

bool StateManager::get_modified_hydration() const {return _modified_hydration;}