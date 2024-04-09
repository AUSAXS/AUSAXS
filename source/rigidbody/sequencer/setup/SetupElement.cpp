/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/setup/SetupElement.h>
#include <rigidbody/sequencer/setup/LoadExistingElement.h>
#include <rigidbody/sequencer/setup/LoadElement.h>
#include <rigidbody/sequencer/setup/AutoConstraintsElement.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/RigidBody.h>

using namespace rigidbody::sequencer;

SetupElement::SetupElement(observer_ptr<Sequencer> owner) : LoopElementCallback(owner) {}

SetupElement& SetupElement::load(const std::vector<std::string>& paths, const std::vector<std::string>& names) {
    owner->_get_elements().push_back(std::make_unique<LoadElement>(static_cast<Sequencer*>(owner), paths, names));
    return *this;
}

SetupElement& SetupElement::load_existing(observer_ptr<RigidBody> rigidbody) {
    owner->_get_elements().push_back(std::make_unique<LoadExistingElement>(static_cast<Sequencer*>(owner), rigidbody));
    return *this;
}

std::unordered_map<std::string, unsigned int>& SetupElement::_get_body_names() {
    return body_names;
}

void SetupElement::_set_active_body(observer_ptr<RigidBody> body) {
    active_body = body;
}

SetupElement& SetupElement::distance_constraint() {
    return *this;
}

SetupElement& SetupElement::fixed_constraint() {
    return *this;
}

SetupElement& SetupElement::generate_linear_constraints() {
    owner->_get_elements().push_back(std::make_unique<AutoConstraintsElement>(static_cast<Sequencer*>(owner), settings::rigidbody::ConstraintGenerationStrategyChoice::Linear));
    return *this;
}

SetupElement& SetupElement::generate_volumetric_constraints() {
    owner->_get_elements().push_back(std::make_unique<AutoConstraintsElement>(static_cast<Sequencer*>(owner), settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric));
    return *this;
}