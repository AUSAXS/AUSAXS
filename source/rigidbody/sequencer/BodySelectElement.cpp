// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::sequencer;

BodySelectElement::BodySelectElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
BodySelectElement::~BodySelectElement() = default;

void BodySelectElement::run() {
    owner->_get_rigidbody()->set_body_select_manager(strategy);
}