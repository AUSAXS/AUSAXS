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
#include <rigidbody/sequencer/OptimizeStepElement.h>

#include <iostream>

using namespace rigidbody::sequencer;

LoopElement::LoopElement() {}

LoopElement::LoopElement(LoopElement* owner) : owner(owner) {}

LoopElement::LoopElement(LoopElement* owner, unsigned int repeats) : owner(owner), iterations(repeats) {}

LoopElement::~LoopElement() = default;

std::shared_ptr<fitter::Fit> LoopElement::execute() {
    std::cout << "LoopElement::execute()" << std::endl;
    return owner->execute(); // propagate upwards to the main Sequencer
}

LoopElement& LoopElement::loop(unsigned int repeats) {
    std::cout << "LoopElement::loop(" << repeats << ")" << std::endl;
    LoopElement& ref = *inner_loops.emplace_back(std::make_unique<LoopElement>(owner, repeats));
    return ref;
}

ParameterElement& LoopElement::parameter_strategy(std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy) {
    std::cout << "LoopElement::parameter_strategy()" << std::endl;
    elements.push_back(std::make_unique<ParameterElement>(this, std::move(strategy)));
    return *static_cast<ParameterElement*>(elements.back().get());
}

BodySelectElement& LoopElement::body_select_strategy(std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) {
    std::cout << "LoopElement::body_select_strategy()" << std::endl;
    elements.push_back(std::make_unique<BodySelectElement>(this, std::move(strategy)));
    return *static_cast<BodySelectElement*>(elements.back().get());
}

TransformElement& LoopElement::transform_strategy(std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) {
    std::cout << "LoopElement::transform_strategy()" << std::endl;
    elements.push_back(std::make_unique<TransformElement>(this, std::move(strategy)));
    return *static_cast<TransformElement*>(elements.back().get());
}

void LoopElement::run() {
    std::cout << "LoopElement::run()" << std::endl;

    for (unsigned int i = 0; i < iterations; i++) {
        for (auto& element : elements) {
            element->run();
        }
    }
}

LoopElement& LoopElement::optimize() {
    std::cout << "LoopElement::optimize()" << std::endl;
    elements.push_back(std::make_unique<OptimizeStepElement>(this));
    return *this;
}

LoopElement& LoopElement::end() {
    std::cout << "LoopElement::end()" << std::endl;
    return *owner;
}