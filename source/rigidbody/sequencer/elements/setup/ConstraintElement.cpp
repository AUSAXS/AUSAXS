// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/ConstraintElement.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::sequencer;

ConstraintElement::ConstraintElement(observer_ptr<Sequencer> owner, std::unique_ptr<rigidbody::constraints::Constraint> constraint) {
    owner->_get_rigidbody()->get_constraint_manager()->add_constraint(std::move(constraint));
}

ConstraintElement::ConstraintElement(observer_ptr<Sequencer> owner, const std::string& body1, const std::string& body2, bool center_mass) {
    owner->_get_rigidbody()->get_constraint_manager()->add_constraint(
        std::make_unique<rigidbody::constraints::DistanceConstraint>(
            owner->_get_rigidbody(),
            owner->setup()->_get_body_names().at(body1),
            owner->setup()->_get_body_names().at(body2),
            center_mass
        )
    );
}

ConstraintElement::ConstraintElement(observer_ptr<Sequencer> owner, const std::string& body1, const std::string& body2, unsigned int iatom1, unsigned int iatom2) {
    owner->_get_rigidbody()->get_constraint_manager()->add_constraint(
        std::make_unique<rigidbody::constraints::DistanceConstraint>(
            owner->_get_rigidbody(),
            owner->setup()->_get_body_names().at(body1),
            owner->setup()->_get_body_names().at(body2),
            iatom1,
            iatom2
        )
    );
}

void ConstraintElement::run() {}