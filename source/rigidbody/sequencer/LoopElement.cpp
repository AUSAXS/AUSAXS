/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/ParameterElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/TransformElement.h>
#include <rigidbody/sequencer/OptimizeStepElement.h>
#include <rigidbody/sequencer/EveryNStepElement.h>
#include <rigidbody/sequencer/SaveElement.h>

using namespace rigidbody::sequencer;

LoopElement::LoopElement(observer_ptr<LoopElement> owner, unsigned int repeats) : iterations(repeats), owner(owner) {
    total_loop_count *= repeats;
}

LoopElement::~LoopElement() = default;

std::shared_ptr<fitter::Fit> LoopElement::execute() {
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