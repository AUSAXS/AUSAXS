/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/setup/SetupElement.h>
#include <rigidbody/sequencer/setup/LoadExistingElement.h>
#include <rigidbody/sequencer/setup/LoadElement.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody::sequencer;

SetupElement::SetupElement(observer_ptr<Sequencer> owner) : LoopElementCallback(owner) {}

SetupElement& SetupElement::load(const std::vector<std::string>& paths, const std::vector<std::string>& names) {
    owner->_get_elements().push_back(std::make_unique<LoadElement>(static_cast<Sequencer*>(owner), paths));
    body_names = names;
    return *this;
}

SetupElement& SetupElement::load_existing(observer_ptr<RigidBody> rigidbody) {
    owner->_get_elements().push_back(std::make_unique<LoadExistingElement>(static_cast<Sequencer*>(owner), rigidbody));
    return *this;
}