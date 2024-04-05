/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>

#include <iostream>

using namespace rigidbody::sequencer;

ParameterElement::ParameterElement(LoopElement* owner, std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}

ParameterElement::~ParameterElement() = default;

void ParameterElement::run() {
    std::cout << "ParameterElement::apply()" << std::endl;
}

ParameterElement& ParameterElement::decay_strategy(std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> strategy) {
    std::cout << "ParameterElement::decay_strategy()" << std::endl;
    this->strategy->set_decay_strategy(std::move(strategy));
    return *this;
}

ParameterElement& ParameterElement::max_rotation_angle(double radians) {
    std::cout << "ParameterElement::max_rotation_angle()" << std::endl;
    strategy->set_max_rotation_angle(radians);
    return *this;
}

ParameterElement& ParameterElement::max_translation_distance(double distance) {
    std::cout << "ParameterElement::max_translation_distance()" << std::endl;
    strategy->set_max_translation_distance(distance);
    return *this;
}