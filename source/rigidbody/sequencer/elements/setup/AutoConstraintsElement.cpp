// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/AutoConstraintsElement.h>
#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <settings/RigidBodySettings.h>

using namespace ausaxs::rigidbody::sequencer;

AutoConstraintsElement::AutoConstraintsElement(observer_ptr<Sequencer> owner, settings::rigidbody::ConstraintGenerationStrategyChoice strategy) : owner(owner), strategy(strategy) {}

void AutoConstraintsElement::run() {
    if (owner->_get_rigidbody() == nullptr) {throw std::runtime_error("AutoConstraintsElement::run: No body is currently loaded.");}
    owner->_get_rigidbody()->constraints->generate_constraints(rigidbody::factory::generate_constraints(owner->_get_rigidbody()->constraints.get(), strategy));
}