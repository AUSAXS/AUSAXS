// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/TransformElement.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::sequencer;

TransformElement::TransformElement(observer_ptr<LoopElement> owner, std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) : LoopElementCallback(owner), strategy(std::move(strategy)) {}
TransformElement::~TransformElement() = default;

void TransformElement::run() {
    owner->_get_rigidbody()->set_transform_manager(strategy);
}