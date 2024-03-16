/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/DecayElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

ParameterElement::ParameterElement(LoopElement* owner) : LoopElementCallback(owner), owner(owner), strategy(settings::rigidbody::parameter_generation_strategy) {
    initialize();
}

ParameterElement::ParameterElement(LoopElement* owner, settings::rigidbody::ParameterGenerationStrategyChoice strategy) : LoopElementCallback(owner), owner(owner), strategy(strategy) {
    initialize();
}

ParameterElement::~ParameterElement() = default;

void ParameterElement::apply() {
    std::cout << "ParameterElement::apply()" << std::endl;
}

DecayElement& ParameterElement::decay_strategy(const settings::rigidbody::DecayStrategyChoice& strategy) {
    std::cout << "ParameterElement::decay_strategy()" << std::endl;
    decay_element = std::make_unique<DecayElement>(this, strategy);
    return *decay_element;
}

ParameterElement& ParameterElement::amplitude(double amplitude) {
    std::cout << "ParameterElement::amplitude()" << std::endl;
    return *this;
}

void ParameterElement::initialize() {
    decay_element = std::make_unique<DecayElement>(this);
}