#include <data/StateManager.h>

#include <memory>
#include <iostream>

StateManager::Signaller::Signaller(unsigned int id, StateManager* const owner) : owner(owner), id(id) {}

void StateManager::Signaller::state_change() const {
    owner->modified(id);
    std::cout << "Signalling object is bound, and a signal was sent." << std::endl;
}

StateManager::UnboundSignaller::UnboundSignaller() : Signaller(0, nullptr) {}

void StateManager::UnboundSignaller::state_change() const {std::cout << "Signalling object is not bound to anything." << std::endl;}

StateManager::StateManager(unsigned int size) : size(size), _modified(size, true), _modified_hydration(true) {
    for (unsigned int i = 0; i < size; i++) {
        probes.push_back(std::make_shared<Signaller>(i, this));
    }
}

void StateManager::modified_all() {
    _modified = std::vector<bool>(size, true);
}

void StateManager::modified(const int i) {
    _modified[i] = true;
}

void StateManager::modified_hydration_layer() {
    _modified_hydration = true;
}

void StateManager::reset() {
    _modified = std::vector<bool>(size, false);
    _modified_hydration = false;
}

std::shared_ptr<StateManager::Signaller> StateManager::get_probe(unsigned int i) {return probes[i];}

std::vector<bool> StateManager::get_modified_bodies() const {return _modified;}

bool StateManager::is_modified(unsigned int i) {return _modified[i];}

bool StateManager::get_modified_hydration() const {return _modified_hydration;}