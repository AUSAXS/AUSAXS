#include <rigidbody/sequencer/TransformElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

TransformElement::TransformElement(LoopElement* owner) : LoopElementCallback(owner), strategy(settings::rigidbody::transform_strategy) {}
TransformElement::TransformElement(LoopElement* owner, settings::rigidbody::TransformationStrategyChoice strategy) : LoopElementCallback(owner), strategy(strategy) {}
TransformElement::~TransformElement() = default;

void TransformElement::apply() {
    std::cout << "TransformElement::apply()" << std::endl;
}