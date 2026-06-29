// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/generation/BackboneConstraints.h>
#include <rigidbody/constraints/generation/NoConstraints.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace rigidbody::constraints;

std::unique_ptr<ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(observer_ptr<const constraints::ConstraintManager> manager) {
    return generate_constraints(manager, settings::rigidbody::constraint_generation_strategy);
}

std::unique_ptr<ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(
    observer_ptr<const constraints::ConstraintManager> manager, const settings::rigidbody::ConstraintGenerationStrategyChoice& choice
) {
    switch (choice) {
        case settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone:
            return std::make_unique<BackboneConstraints>(manager);
        case settings::rigidbody::ConstraintGenerationStrategyChoice::None:
            return std::make_unique<NoConstraints>(manager);
        default: 
            throw except::unexpected("rigidbody::factory::generate_constraints: Unknown constraint generation strategy choice. Did you forget to add it to the switch statement?");
    }
}