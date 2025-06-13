// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/LoadExistingElement.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::sequencer;

LoadExistingElement::LoadExistingElement(observer_ptr<Sequencer> owner, observer_ptr<RigidBody> rigidbody) : owner(owner), rigidbody(rigidbody) {
    owner->setup()->_set_active_body(rigidbody);
}

void LoadExistingElement::run() {
    owner->_get_rigidbody() = rigidbody;
}