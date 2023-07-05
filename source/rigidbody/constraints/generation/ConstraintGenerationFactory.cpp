#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/generation/VolumetricConstraints.h>
#include <rigidbody/constraints/generation/NoConstraints.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

std::unique_ptr<rigidbody::ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(const ConstraintManager* manager) {
    return generate_constraints(manager, settings::rigidbody::constraint_generation_strategy);
}

std::unique_ptr<rigidbody::ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(const ConstraintManager* manager, const settings::rigidbody::ConstraintGenerationStrategyChoice& choice) {
    switch (choice) {
        case settings::rigidbody::ConstraintGenerationStrategyChoice::Linear:
            return std::make_unique<rigidbody::LinearConstraints>(manager);
        case settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric:
            return std::make_unique<rigidbody::VolumetricConstraints>(manager);
        case settings::rigidbody::ConstraintGenerationStrategyChoice::None: 
            return std::make_unique<rigidbody::NoConstraints>(manager);
        default: 
            throw except::unexpected("rigidbody::factory::generate_constraints: Unknown constraint generation strategy choice. Did you forget to add it to the switch statement?");
    }
}