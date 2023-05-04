#include <rigidbody/sequencer/ParameterElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

ParameterElement::ParameterElement(LoopElement* owner) 
    : owner(owner), LoopElementCallback(owner), strategy(settings::rigidbody::parameter_generation_strategy) {}

ParameterElement::ParameterElement(LoopElement* owner, settings::rigidbody::ParameterGenerationStrategyChoice strategy) 
    : owner(owner), LoopElementCallback(owner), strategy(strategy) {}

ParameterElement::~ParameterElement() = default;

void ParameterElement::apply() {
    std::cout << "ParameterElement::apply()" << std::endl;
}

DecayElement& ParameterElement::decay_strategy(settings::rigidbody::DecayStrategyChoice strategy) {
    std::cout << "ParameterElement::decay_strategy()" << std::endl;
    decay_element = std::make_unique<DecayElement>(owner, strategy);
    return *decay_element;
}

ParameterElement& ParameterElement::amplitude(double amplitude) {
    std::cout << "ParameterElement::amplitude()" << std::endl;
    return *this;
}
