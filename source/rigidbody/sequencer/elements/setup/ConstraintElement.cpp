// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/ConstraintElement.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::sequencer;

ConstraintElement::ConstraintElement(observer_ptr<Sequencer> owner, std::unique_ptr<rigidbody::constraints::Constraint> constraint) {
    owner->_get_rigidbody()->constraints->add_constraint(std::move(constraint));
}

void ConstraintElement::run() {}