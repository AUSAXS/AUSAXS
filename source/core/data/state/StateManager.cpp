/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>

#ifdef DEBUG
    #include <stdexcept>
    #include <string>
#endif

using namespace ausaxs;
using namespace ausaxs::state;

StateManager::StateManager(unsigned int size) : _size(size), _externally_modified(size, true), _internally_modified(size, true), _modified_hydration(true) {
    for (unsigned int i = 0; i < size; i++) {
        probes.push_back(std::make_shared<signaller::BoundSignaller>(i, this));
    }
}

void StateManager::externally_modified_all() {
    _externally_modified = std::vector<bool>(size(), true);
}

void StateManager::internally_modified_all() {
    _internally_modified = std::vector<bool>(size(), true);
}

void StateManager::externally_modified(unsigned int i) {
    #if DEBUG
        if (i >= size()) {
            throw std::out_of_range("StateManager::externally_modified: index out of range (" + std::to_string(i) + " >= " + std::to_string(size()) + ")");
        }
    #endif
    _externally_modified[i] = true;
}

void StateManager::internally_modified(unsigned int i) {
    #if DEBUG
        if (i >= size()) {
            throw std::out_of_range("StateManager::internally_modified: index out of range (" + std::to_string(i) + " >= " + std::to_string(size()) + ")");
        }
    #endif
    _internally_modified[i] = true;
}

void StateManager::modified_hydration_layer() {
    _modified_hydration = true;
}

void StateManager::reset_to_false() {
    _internally_modified = std::vector<bool>(size(), false);
    _externally_modified = std::vector<bool>(size(), false);
    _modified_hydration = false;
}

void StateManager::set_probe(unsigned int i, std::shared_ptr<signaller::Signaller> probe) {
    #if DEBUG
        if (probes.size() <= i) {
            throw std::out_of_range("StateManager::set_probe: index out of range (" + std::to_string(i) + " >= " + std::to_string(probes.size()) + ")");
        }
    #endif
    probes[i] = std::move(probe);
}

std::shared_ptr<signaller::Signaller> StateManager::get_probe(unsigned int i) {
    #if DEBUG
        if (probes.size() <= i) {
            throw std::out_of_range("StateManager::get_probe: index out of range (" + std::to_string(i) + " >= " + std::to_string(probes.size()) + ")");
        }
    #endif
    return probes[i];
}

std::vector<std::shared_ptr<signaller::Signaller>> StateManager::get_probes() {return probes;}

const std::vector<bool>& StateManager::get_externally_modified_bodies() const {return _externally_modified;}

const std::vector<bool>& StateManager::get_internally_modified_bodies() const {return _internally_modified;}

bool StateManager::is_externally_modified(unsigned int i) {return _externally_modified[i];}

bool StateManager::is_internally_modified(unsigned int i) {return _internally_modified[i];}

bool StateManager::get_modified_hydration() const {return _modified_hydration;}

std::size_t StateManager::size() const {return _size;}