/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::sequencer;

ParameterElement::ParameterElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}

ParameterElement::~ParameterElement() = default;

void ParameterElement::run() {
    owner->_get_rigidbody()->set_parameter_manager(strategy);
}

ParameterElement& ParameterElement::decay_strategy(std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> strategy) {
    this->strategy->set_decay_strategy(std::move(strategy));
    return *this;
}

ParameterElement& ParameterElement::max_rotation_angle(double radians) {
    strategy->set_max_rotation_angle(radians);
    return *this;
}

ParameterElement& ParameterElement::max_translation_distance(double distance) {
    strategy->set_max_translation_distance(distance);
    return *this;
}