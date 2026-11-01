// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/EveryNStepElement.h>

#include <rigidbody/sequencer/detail/parse_error.h>
#include <utility/StringUtils.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

EveryNStepElement::EveryNStepElement(observer_ptr<LoopElement> owner, int n) : LoopElement(owner, 1), n(n) {
    assert(0 < n && "EveryNStepElement: n must be positive; run() and the step count both divide by it.");
}

EveryNStepElement::~EveryNStepElement() = default;

void EveryNStepElement::run() {
    if (++loop_counter % n == 0) {
        for (auto& e : elements) {
            e->run();
        }
    }
}

int EveryNStepElement::_get_step_size() const {
    return n;
}

std::vector<std::string> EveryNStepElement::_valid_arguments() {
    return {};
}

InlineSignature EveryNStepElement::_valid_inline_arguments() {
    return {.names = {"steps"}, .min = 1, .max = 1};
}

// every [steps] - opens a block run once every [steps] iterations
std::unique_ptr<GenericElement> EveryNStepElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) { // NOLINT
    if (!utility::isinteger(args.inlined[0])) {throw except::parse_error("every_n_step", "Expected an integer value for the number of steps, but got \"" + args.inlined[0] + "\".");}
    int n = search::parse_int(args.inlined[0]);
    if (n <= 0) {throw except::parse_error("every_n_step", "The number of steps must be positive, but got " + args.inlined[0] + ".");} // run() takes it modulo
    return std::make_unique<EveryNStepElement>(owner, n);
}