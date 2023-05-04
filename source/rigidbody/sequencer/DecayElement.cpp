#include <rigidbody/sequencer/DecayElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

DecayElement::DecayElement(ParameterElement* owner) 
    : ParameterElementCallback(owner), strategy(settings::rigidbody::decay_strategy) {}
DecayElement::DecayElement(ParameterElement* owner, settings::rigidbody::DecayStrategyChoice strategy) 
    : ParameterElementCallback(owner), strategy(strategy) {}
DecayElement::~DecayElement() = default;

void DecayElement::apply() {
    std::cout << "DecayElement::apply()" << std::endl;
}