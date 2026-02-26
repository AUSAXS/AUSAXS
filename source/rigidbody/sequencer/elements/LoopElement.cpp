// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/CopyLoopElement.h>
#include <rigidbody/sequencer/elements/ParameterElement.h>
#include <rigidbody/sequencer/elements/BodySelectElement.h>
#include <rigidbody/sequencer/elements/TransformElement.h>
#include <rigidbody/sequencer/elements/OptimizeStepElement.h>
#include <rigidbody/sequencer/elements/EveryNStepElement.h>
#include <rigidbody/sequencer/elements/SaveElement.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>

#include <cassert>

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

observer_ptr<rigidbody::Rigidbody> LoopElement::_get_rigidbody() const {
    assert(owner != nullptr && "LoopElement::_get_rigidbody: Owner is null.");
    return owner->_get_rigidbody();
}

observer_ptr<data::Molecule> LoopElement::_get_molecule() const {
    assert(owner != nullptr && "LoopElement::_get_molecule: Owner is null.");
    return owner->_get_molecule();
}

observer_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> LoopElement::_get_best_conf() const {
    assert(owner != nullptr && "LoopElement::_get_best_conf: Owner is null.");
    return owner->_get_best_conf();
}

observer_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> LoopElement::_get_current_conf() const {
    return &_get_rigidbody()->conformation->absolute_parameters;
}

observer_ptr<LoopElement> LoopElement::_get_owner() const {
    assert(owner != nullptr && "LoopElement::_get_owner: Owner is null.");
    return owner;
}

observer_ptr<const Sequencer> LoopElement::_get_sequencer() const {
    assert(owner != nullptr && "LoopElement::_get_sequencer: Owner is null.");
    return owner->_get_sequencer();
}

observer_ptr<Sequencer> LoopElement::_get_sequencer() {
    assert(owner != nullptr && "LoopElement::_get_sequencer: Owner is null.");
    return owner->_get_sequencer();
}

std::vector<std::unique_ptr<GenericElement>>& LoopElement::_get_elements() {
    return elements;
}

unsigned int LoopElement::_get_loop_iterations() const {
    return iterations;
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

unsigned int LoopElement::_get_current_iteration() {
    return global_counter;
}

unsigned int LoopElement::_get_total_iterations() {
    return total_loop_count;
}

void LoopElement::_add_total_iterations(unsigned int n) {
    total_loop_count += n;
}

void LoopElement::_parse::validate(observer_ptr<LoopElement> owner, const ParsedArgs& args) {}

std::unique_ptr<GenericElement> LoopElement::_parse::create(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    static std::unordered_map<std::string, observer_ptr<LoopElement>> loop_names;
    static observer_ptr<LoopElement> last_loop_element = nullptr;

    auto deduce_iteration_count = [&]() -> int {
        // find the last parameter element by traversing backwards and upwards through the owner chain
        auto find_last_parameter_element = [&]() -> observer_ptr<ParameterElement> {
            observer_ptr<LoopElement> current = owner;
            int escape_counter = 0;
            while (current != nullptr) {
                for (auto& e : current->_get_elements()) {
                    if (auto* parameter_element = dynamic_cast<ParameterElement*>(e.get())) {
                        return parameter_element;
                    }
                }
                current = current->_get_owner();
                if (100 < ++escape_counter) {throw std::runtime_error("LoopElement::_parse::create: owner chain too long while searching for last parameter element.");}
            }
            return nullptr;
        };

        auto* last_parameter_element = find_last_parameter_element();
        if (!last_parameter_element) {throw except::parse_error("loop", "Could not deduce number of iterations.");}
        int iterations = last_parameter_element->get_parameter_strategy()->get_decay_strategy()->get_iterations();
        return iterations;
    };

    if (!args.named.empty()) {throw except::parse_error("loop", "Unexpected named argument.");}
    if (args.inlined.values.empty()) { // no args - try to deduce iteration count
        return std::make_unique<LoopElement>(owner, deduce_iteration_count());
    } else if (args.inlined.values.size() == 1) {
        // five options: 1) "[iteration]" 2) "[name]" 3) "[name] [iteration]" 4) "[duplicate]" 5) "[duplicate] [name]"
        if (args.inlined.values.size() == 1) { // option 1, 2, 4
            try { // check option 1
                int iterations = std::stoi(args.inlined.values[0]);
                return std::make_unique<LoopElement>(owner, iterations);
            } catch (std::exception&) {
                const auto& name = args.inlined.values[0];

                // check option 4
                if (name == "duplicate" || name == "copy") {
                    if (loop_names.empty()) {throw except::parse_error("loop", args.inlined, "No previous loop found to duplicate.");}
                    return std::make_unique<CopyLoopElement>(owner, last_loop_element);
                }

                // else it must be option 2
                if (loop_names.contains(name)) {throw except::parse_error("loop", args.inlined, "Loop name \"" + name + "\" already exists.");}
                auto loop = std::make_unique<LoopElement>(owner, deduce_iteration_count());
                loop_names[name] = loop.get();
                last_loop_element = loop.get();
                return loop;
            }
        } else if (args.inlined.values.size() == 2) { // option 3, 5
            // check option 5
            if (args.inlined.values[0] == "duplicate" || args.inlined.values[0] == "copy") {
                const auto& name = args.inlined.values[1];
                if (!loop_names.contains(name)) {throw except::parse_error("loop", args.inlined, "Target loop name \"" + name + "\" does not exist.");}
                return std::make_unique<CopyLoopElement>(owner, loop_names.at(name));
            }

            // else it must be option 3
            try {
                int iterations = std::stoi(args.inlined.values[1]);
                const auto& name = args.inlined.values[0];
                if (loop_names.contains(name)) {throw except::parse_error("loop", args.inlined, "Loop name \"" + name + "\" already exists.");}
                auto loop = std::make_unique<LoopElement>(owner, iterations);
                loop_names[name] = loop.get();
                last_loop_element = loop.get();
                return loop;
            } catch (std::exception&) {
                throw except::parse_error("loop", args.inlined, "Could not determine deduce number of iterations.");
            }
        }
    }
    throw except::parse_error("loop", args.inlined, "Invalid arguments.");
}