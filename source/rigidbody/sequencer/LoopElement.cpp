#include <rigidbody/sequencer/LoopElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

LoopElement::LoopElement(LoopElement* owner) : owner(owner) {}
LoopElement::LoopElement(unsigned int iterations) : iterations(iterations) {}
                
void LoopElement::execute() {
    std::cout << "LoopElement::execute()" << std::endl;
    parameter_element->apply();
    body_select_element->apply();
    transform_element->apply();
    for (unsigned int i = 0; i < iterations; i++) {
        this->run();
        for (auto& loop : inner_loops) {loop->execute();}
    }
}

LoopElement& LoopElement::loop(unsigned int repeats) {
    std::cout << "LoopElement::loop(" << repeats << ")" << std::endl;
    return *inner_loops.emplace_back(std::make_unique<LoopElement>(repeats));
}

ParameterElement& LoopElement::parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice strategy) {
    std::cout << "LoopElement::parameter_strategy()" << std::endl;
    parameter_element = std::make_unique<ParameterElement>(this, strategy);
    return *parameter_element;
}

BodySelectElement& LoopElement::body_select_strategy(settings::rigidbody::BodySelectStrategyChoice strategy) {
    std::cout << "LoopElement::body_select_strategy()" << std::endl;
    body_select_element = std::make_unique<BodySelectElement>(this, strategy);
    return *body_select_element;
}

TransformElement& LoopElement::transform_strategy(settings::rigidbody::TransformationStrategyChoice strategy) {
    std::cout << "LoopElement::transform_strategy()" << std::endl;
    transform_element = std::make_unique<TransformElement>(this, strategy);
    return *transform_element;
}

void LoopElement::run() {
    std::cout << "LoopElement::run()" << std::endl;
}