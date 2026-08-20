// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/CopyLoopElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>

using namespace ausaxs::rigidbody::sequencer; 

CopyLoopElement::CopyLoopElement(observer_ptr<LoopElement> owner, observer_ptr<LoopElement> target) : owner(owner), target(target) {}

CopyLoopElement::~CopyLoopElement() = default;

void CopyLoopElement::run() {
    target->run();
}

ausaxs::observer_ptr<LoopElement> CopyLoopElement::_get_target() const {
    return target;
}