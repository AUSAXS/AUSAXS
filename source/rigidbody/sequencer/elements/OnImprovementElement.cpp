// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/OnImprovementElement.h>
#include <rigidbody/sequencer/elements/OptimizeStepElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>

using namespace ausaxs::rigidbody::sequencer;

OnImprovementElement::OnImprovementElement(observer_ptr<OptimizeStepElement> owner)
    : LoopElement(owner, 1), owner(owner) {}

OnImprovementElement::~OnImprovementElement() = default;

void OnImprovementElement::run() {
    // Only run children if the parent optimization step was accepted
    if (owner->was_accepted()) {
        for (auto& e : elements) {
            e->run();
        }
    }
}

std::vector<std::string> OnImprovementElement::_valid_arguments() {
    return {};
}

InlineSignature OnImprovementElement::_valid_inline_arguments() {
    return {.names = {}, .min = 0, .max = 0};
}

std::unique_ptr<GenericElement> OnImprovementElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&&) {
    observer_ptr<OptimizeStepElement> optimize_step = nullptr;
    if (optimize_step = dynamic_cast<OptimizeStepElement*>(owner); !optimize_step) {
        throw except::parse_error("on_improvement", "\"on_improvement\" must be inside an \"optimize_step\" block.");
    }
    return std::make_unique<OnImprovementElement>(optimize_step);
}