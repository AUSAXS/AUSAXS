/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/setup/LoadElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <settings/GeneralSettings.h>

#include <iostream>

using namespace rigidbody::sequencer;

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& path, const std::vector<std::string>& body_names) : owner(owner) {
    rigidbody = std::make_unique<RigidBody>(data::Molecule(path));
    if (!body_names.empty() && body_names.size() != rigidbody->body_size()) {throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");}
    for (unsigned int i = 0; i < rigidbody->body_size(); ++i) {
        owner->_get_body_names().emplace(body_names.empty() ? "b" + std::to_string(i) : body_names[i], i);
    }
    owner->_set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->body_size() << " bodies from " << path.size() << " files." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& splits, const std::vector<std::string>& body_names) : owner(owner) {
    rigidbody = std::make_unique<RigidBody>(rigidbody::BodySplitter::split(path, splits));
    if (!body_names.empty() && body_names.size() != rigidbody->body_size()) {throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");}
    for (unsigned int i = 0; i < rigidbody->body_size(); ++i) {
        owner->_get_body_names().emplace(body_names.empty() ? "b" + std::to_string(i) : body_names[i], i);
    }
    owner->_set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->body_size() << " bodies from \"" << path.size() << "\"." << std::endl;
    }
}

void LoadElement::run() {
    owner->_get_rigidbody() = rigidbody.get();
}