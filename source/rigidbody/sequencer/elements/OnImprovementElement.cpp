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

std::unique_ptr<GenericElement> OnImprovementElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("on_improvement", "Unexpected named argument.");}
    if (!args.inlined.empty()) {throw except::parse_error("on_improvement", "Unexpected inline argument.");}

    observer_ptr<OptimizeStepElement> optimize_step = nullptr;
    if (optimize_step = dynamic_cast<OptimizeStepElement*>(owner); optimize_step) {
        throw except::parse_error("on_improvement", "\"on_improvement\" must be inside an \"optimize_step\" block.");
    }
    return std::make_unique<OnImprovementElement>(optimize_step);
}