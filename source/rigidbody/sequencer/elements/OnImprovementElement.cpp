// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/OnImprovementElement.h>
#include <rigidbody/sequencer/elements/OptimizeStepElement.h>
#include <rigidbody/detail/Configuration.h>

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