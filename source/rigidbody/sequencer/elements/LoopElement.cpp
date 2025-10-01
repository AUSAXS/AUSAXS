// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/elements/ParameterElement.h>
#include <rigidbody/sequencer/elements/BodySelectElement.h>
#include <rigidbody/sequencer/elements/TransformElement.h>
#include <rigidbody/sequencer/elements/OptimizeStepElement.h>
#include <rigidbody/sequencer/elements/EveryNStepElement.h>
#include <rigidbody/sequencer/elements/SaveElement.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

LoopElement::LoopElement(observer_ptr<LoopElement> owner, unsigned int repeats) : iterations(repeats), owner(owner) {
    if (iterations == 1) {return;}
    int this_will_run = iterations;
    auto next_owner = _get_owner();
    int escape_counter = 0;
    while (dynamic_cast<Sequencer*>(next_owner) == nullptr) {
        if (100 < ++escape_counter) {throw std::runtime_error("LoopElement::LoopElement: owner chain too long");}
        this_will_run *= next_owner->iterations;
        next_owner = next_owner->_get_owner();
    }
    total_loop_count += this_will_run;
}

LoopElement::~LoopElement() = default;

std::shared_ptr<fitter::FitResult> LoopElement::execute() {
    return owner->execute(); // propagate upwards to the main Sequencer
}

LoopElement& LoopElement::loop(unsigned int repeats) {
    elements.push_back(std::make_unique<LoopElement>(this, repeats));
    return *static_cast<LoopElement*>(elements.back().get());
}

ParameterElement& LoopElement::parameter_strategy(std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy) {
    elements.push_back(std::make_unique<ParameterElement>(this, std::move(strategy)));
    return *static_cast<ParameterElement*>(elements.back().get());
}

BodySelectElement& LoopElement::body_select_strategy(std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) {
    elements.push_back(std::make_unique<BodySelectElement>(this, std::move(strategy)));
    return *static_cast<BodySelectElement*>(elements.back().get());
}

TransformElement& LoopElement::transform_strategy(std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) {
    elements.push_back(std::make_unique<TransformElement>(this, std::move(strategy)));
    return *static_cast<TransformElement*>(elements.back().get());
}

void LoopElement::run() {
    for (unsigned int i = 0; i < iterations; ++i) {
        ++global_counter;
        for (auto& element : elements) {
            element->run();
        }
    }
}

observer_ptr<rigidbody::RigidBody> LoopElement::_get_rigidbody() const {
    return owner->_get_rigidbody();
}

observer_ptr<rigidbody::detail::BestConf> LoopElement::_get_best_conf() const {
    return owner->_get_best_conf();
}

observer_ptr<LoopElement> LoopElement::_get_owner() const {
    return owner;
}

observer_ptr<const Sequencer> LoopElement::_get_sequencer() const {
    return owner->_get_sequencer();
}

std::vector<std::unique_ptr<GenericElement>>& LoopElement::_get_elements() {
    return elements;
}

OptimizeStepElement& LoopElement::optimize() {
    elements.push_back(std::make_unique<OptimizeStepElement>(this));
    return *static_cast<OptimizeStepElement*>(elements.back().get());
}

LoopElement& LoopElement::end() {
    return *owner;
}

LoopElement& LoopElement::save(const io::File& path) {
    elements.push_back(std::make_unique<SaveElement>(this, path));
    return *this;
}

EveryNStepElement& LoopElement::every(unsigned int n) {
    elements.push_back(std::make_unique<EveryNStepElement>(this, n));
    return *static_cast<EveryNStepElement*>(elements.back().get());
}