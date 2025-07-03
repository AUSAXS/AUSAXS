// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/setup/AutoConstraintsElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/RigidBody.h>
#include <settings/RigidBodySettings.h>

using namespace ausaxs::rigidbody::sequencer;

AutoConstraintsElement::AutoConstraintsElement(observer_ptr<Sequencer> owner, settings::rigidbody::ConstraintGenerationStrategyChoice strategy) : owner(owner), strategy(strategy) {}

void AutoConstraintsElement::run() {
    if (owner->_get_rigidbody() == nullptr) {throw std::runtime_error("AutoConstraintsElement::run: No body is currently loaded.");}
    owner->_get_rigidbody()->get_constraint_manager()->generate_constraints(rigidbody::factory::generate_constraints(owner->_get_rigidbody()->get_constraint_manager().get(), strategy));
}