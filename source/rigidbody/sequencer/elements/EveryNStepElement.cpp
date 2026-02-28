// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/EveryNStepElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <utility/StringUtils.h>

using namespace ausaxs::rigidbody::sequencer;

EveryNStepElement::EveryNStepElement(observer_ptr<LoopElement> owner, unsigned int n) : LoopElement(owner, 1), n(n), loop_counter(0) {}

EveryNStepElement::~EveryNStepElement() = default;

void EveryNStepElement::run() {
    if (++loop_counter % n == 0) {
        for (auto& e : elements) {
            e->run();
        }
    }
}

std::unique_ptr<GenericElement> EveryNStepElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("every_n_step", "Unexpected named argument.");}
    if (args.inlined.size() != 1) {throw except::parse_error("every_n_step", "Expected exactly one inline argument.");}
    if (!utility::isinteger(args.inlined[0])) {throw except::parse_error("every_n_step", "Expected an integer value for the number of steps, but got \"" + args.inlined[0] + "\".");}
    return std::make_unique<EveryNStepElement>(owner, std::stoi(args.inlined[0]));
}