// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/EveryNStepElement.h>

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