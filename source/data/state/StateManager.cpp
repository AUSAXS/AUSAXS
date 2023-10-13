#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>

#include <iostream>

using namespace state;

StateManager::StateManager(unsigned int size) : size(size), _externally_modified(size, true), _internally_modified(size, true), _modified_hydration(true) {
    for (unsigned int i = 0; i < size; i++) {
        probes.push_back(std::make_shared<signaller::BoundSignaller>(i, this));
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

void StateManager::set_probe(unsigned int i, std::shared_ptr<signaller::Signaller> probe) {probes[i] = probe;}

std::shared_ptr<signaller::Signaller> StateManager::get_probe(unsigned int i) {return probes[i];}

std::vector<std::shared_ptr<signaller::Signaller>> StateManager::get_probes() {return probes;}

const std::vector<bool>& StateManager::get_externally_modified_bodies() const {return _externally_modified;}

const std::vector<bool>& StateManager::get_internally_modified_bodies() const {return _internally_modified;}

bool StateManager::is_externally_modified(unsigned int i) {return _externally_modified[i];}

bool StateManager::is_internally_modified(unsigned int i) {return _internally_modified[i];}

bool StateManager::get_modified_hydration() const {return _modified_hydration;}