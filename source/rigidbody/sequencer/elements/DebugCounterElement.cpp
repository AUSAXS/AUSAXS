// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/DebugCounterElement.h>

#include <rigidbody/sequencer/elements/LoopElement.h>

using namespace ausaxs::rigidbody::sequencer;

DebugCounterElement::DebugCounterElement(std::string name) : label(std::move(name)) {}

DebugCounterElement::~DebugCounterElement() = default;

void DebugCounterElement::run() {
    recorded.push_back(LoopElement::_get_current_iteration());
}

int DebugCounterElement::count() const {
    return static_cast<int>(recorded.size());
}

const std::vector<int>& DebugCounterElement::iterations() const {
    return recorded;
}

const std::string& DebugCounterElement::name() const {
    return label;
}

std::vector<std::string> DebugCounterElement::_valid_arguments() {
    return {};
}

InlineSignature DebugCounterElement::_valid_inline_arguments() {
    return {.names = {"name"}, .min = 0, .max = 1};
}

// debug_counter [name] - records each run; for tests only
std::unique_ptr<GenericElement> DebugCounterElement::_parse(observer_ptr<LoopElement>, ParsedArgs&& args) { // NOLINT
    return std::make_unique<DebugCounterElement>(args.inlined.empty() ? "" : args.inlined[0]);
}
