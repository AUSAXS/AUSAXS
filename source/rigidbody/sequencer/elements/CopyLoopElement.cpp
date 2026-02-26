// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/CopyLoopElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/Sequencer.h>

using namespace ausaxs::rigidbody::sequencer; 

CopyLoopElement::CopyLoopElement(observer_ptr<LoopElement> owner, observer_ptr<LoopElement> target) : owner(owner), target(target) {
    // Replicate the total_loop_count contribution that LoopElement::LoopElement computes,
    // using the target's iteration count and the owner chain of this element's position.
    int this_will_run = static_cast<int>(target->_get_loop_iterations());
    auto next_owner = owner;
    int escape_counter = 0;
    while (dynamic_cast<Sequencer*>(next_owner) == nullptr) {
        if (100 < ++escape_counter) {throw std::runtime_error("CopyLoopElement::CopyLoopElement: owner chain too long");}
        this_will_run *= static_cast<int>(next_owner->_get_loop_iterations());
        next_owner = next_owner->_get_owner();
    }
    LoopElement::_add_total_iterations(this_will_run);
}

CopyLoopElement::~CopyLoopElement() = default;

void CopyLoopElement::run() {
    target->run();
}