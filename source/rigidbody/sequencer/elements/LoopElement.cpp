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

namespace {
    std::unordered_map<std::string, observer_ptr<LoopElement>> loop_names;
    observer_ptr<LoopElement> last_loop_element = nullptr;
}

LoopElement::LoopElement(observer_ptr<LoopElement> owner, unsigned int repeats) : iterations(repeats), owner(owner) {}

void LoopElement::_reset_counters() {
    total_loop_count = 0;
    global_counter = 0;
}

void LoopElement::_reset_named_loops() {
    loop_names.clear();
    last_loop_element = nullptr;
}

LoopElement::~LoopElement() {
    _reset_counters();
    _reset_named_loops();
}

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
        // checked here rather than between the individual elements so a stopped iteration is never left half-finished
        if (_stop_requested()) {return;}
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

void LoopElement::_request_stop() {
    stop_flag.store(true, std::memory_order_relaxed);
}

bool LoopElement::_stop_requested() {
    return stop_flag.load(std::memory_order_relaxed);
}

void LoopElement::_clear_stop_request() {
    stop_flag.store(false, std::memory_order_relaxed);
}

unsigned int LoopElement::_get_current_iteration() {
    return global_counter;
}

unsigned int LoopElement::_get_total_iterations() {
    return total_loop_count;
}

namespace {
    // Number of optimization steps performed by a single full run of the given loop, where multiplier is the
    // number of times the loop body itself is run. Nested loops multiply the counter by their own iteration
    // count as they are entered, so a step deep in the tree is counted once for every time it is reached.
    unsigned int count_optimization_steps(observer_ptr<LoopElement> loop, unsigned int multiplier, int depth = 0) {
        if (100 < ++depth) {throw std::runtime_error("LoopElement::count_optimization_steps: element tree too deep");}

        unsigned int steps = 0;
        for (auto& e : loop->_get_elements()) {
            // a copy loop runs its target in-place, so it contributes exactly what the target would here
            if (auto* copy = dynamic_cast<CopyLoopElement*>(e.get())) {
                auto* target = copy->_get_target();
                steps += count_optimization_steps(target, multiplier*target->_get_loop_iterations(), depth);
                continue;
            }

            auto* nested = dynamic_cast<LoopElement*>(e.get());
            if (nested == nullptr) {continue;} // only loop-like elements can contain optimization steps

            // an optimize element performs one step itself, and may still hold nested blocks below it
            if (dynamic_cast<OptimizeStepElement*>(nested) != nullptr) {steps += multiplier;}

            unsigned int inner_multiplier = multiplier*nested->_get_loop_iterations();

            // an every-n block only runs its contents on every nth iteration of the surrounding loop
            if (auto* every = dynamic_cast<EveryNStepElement*>(nested)) {inner_multiplier /= every->_get_step_size();}

            steps += count_optimization_steps(nested, inner_multiplier, depth);
        }
        return steps;
    }
}

void LoopElement::_recount_total_iterations(observer_ptr<LoopElement> root) {
    total_loop_count = count_optimization_steps(root, root->iterations);
}

std::vector<std::string> LoopElement::_valid_arguments() {
    return {};
}

InlineSignature LoopElement::_valid_inline_arguments() {
    return {.names = {"name", "iterations"}, .min = 0, .max = 2};
}

std::unique_ptr<GenericElement> LoopElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
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
                // the Sequencer is the top of the chain and has no owner to continue to, so the search ends here
                if (dynamic_cast<Sequencer*>(current) != nullptr) {break;}
                current = current->_get_owner();
                if (100 < ++escape_counter) {throw std::runtime_error("LoopElement::_parse::create: owner chain too long while searching for last parameter element.");}
            }
            return nullptr;
        };

        auto* last_parameter_element = find_last_parameter_element();
        if (!last_parameter_element) {
            throw except::parse_error("loop", "Could not deduce the number of iterations: no preceding \"parameter\" element was found. "
                "Either add one, or state the count explicitly as e.g. \"loop 100\".");
        }
        int iterations = last_parameter_element->get_parameter_strategy()->get_decay_strategy()->get_iterations();
        return iterations;
    };

    if (args.inlined.empty()) { // pattern 1: [] - iteration count deduced from the last parameter element
        return std::make_unique<LoopElement>(owner, deduce_iteration_count());
    } else if (args.inlined.size() == 1) {
        try { // pattern 2: [iterations]
            int iterations = std::stoi(args.inlined[0]);
            return std::make_unique<LoopElement>(owner, iterations);
        } catch (std::exception&) {
            const auto& name = args.inlined[0];

            // pattern 3: [duplicate|copy] - repeats the most recently named loop
            if (name == "duplicate" || name == "copy") {
                if (loop_names.empty()) {throw except::parse_error("loop", args.inlined, "No previous loop found to duplicate.");}
                return std::make_unique<CopyLoopElement>(owner, last_loop_element);
            }

            // pattern 4: [name] - iteration count deduced as in pattern 1
            if (loop_names.contains(name)) {throw except::parse_error("loop", args.inlined, "Loop name \"" + name + "\" already exists.");}
            auto loop = std::make_unique<LoopElement>(owner, deduce_iteration_count());
            loop_names[name] = loop.get();
            last_loop_element = loop.get();
            return loop;
        }
    } else {
        // pattern 5: [duplicate|copy] [name] - repeats the named loop
        if (args.inlined[0] == "duplicate" || args.inlined[0] == "copy") {
            const auto& name = args.inlined[1];
            if (!loop_names.contains(name)) {throw except::parse_error("loop", args.inlined, "Target loop name \"" + name + "\" does not exist.");}
            return std::make_unique<CopyLoopElement>(owner, loop_names.at(name));
        }

        // pattern 6: [name] [iterations]
        try {
            int iterations = std::stoi(args.inlined[1]);
            const auto& name = args.inlined[0];
            if (loop_names.contains(name)) {throw except::parse_error("loop", args.inlined, "Loop name \"" + name + "\" already exists.");}
            auto loop = std::make_unique<LoopElement>(owner, iterations);
            loop_names[name] = loop.get();
            last_loop_element = loop.get();
            return loop;
        } catch (std::exception&) {
            throw except::parse_error("loop", args.inlined, "Could not determine number of iterations.");
        }
    }
}