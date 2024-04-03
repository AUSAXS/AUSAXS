/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/RigidBodyManager.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/TransformElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

LoopElement::LoopElement() : 
    parameter_element(std::make_unique<ParameterElement>(this)), 
    body_select_element(std::make_unique<BodySelectElement>(this)), 
    transform_element(std::make_unique<TransformElement>(this)) 
{}

LoopElement::LoopElement(LoopElement* owner) : 
    owner(owner),
    parameter_element(std::make_unique<ParameterElement>(this)), 
    body_select_element(std::make_unique<BodySelectElement>(this)), 
    transform_element(std::make_unique<TransformElement>(this)) 
{}

LoopElement::LoopElement(LoopElement* owner, unsigned int repeats) : 
    owner(owner),
    iterations(repeats),
    parameter_element(std::make_unique<ParameterElement>(this)), 
    body_select_element(std::make_unique<BodySelectElement>(this)), 
    transform_element(std::make_unique<TransformElement>(this)) 
{}

LoopElement::~LoopElement() = default;

void LoopElement::execute() {
    std::cout << "LoopElement::execute()" << std::endl;
    parameter_element->apply();
    body_select_element->apply();
    transform_element->apply();

    // if this is a root node, run it the given number of times. 
    if (inner_loops.empty()) {
        for (unsigned int i = 0; i < iterations; i++) {
            this->run();
        }
    } else { // else run the inner loops the given number of times.
        for (unsigned int i = 0; i < iterations; i++) {
            for (auto& loop : inner_loops) {loop->execute();}
        }
    }
}

LoopElement& LoopElement::loop(unsigned int repeats) {
    std::cout << "LoopElement::loop(" << repeats << ")" << std::endl;
    LoopElement& ref = *inner_loops.emplace_back(std::make_unique<LoopElement>(owner, repeats));
    ref.body_select_element->strategy = body_select_element->strategy;
    ref.parameter_element->strategy = parameter_element->strategy;
    ref.transform_element->strategy = transform_element->strategy;
    return ref;
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
    rigidbody->optimize_step();
}

LoopElement& LoopElement::end() {
    std::cout << "LoopElement::end()" << std::endl;
    return *owner;
}