#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/generation/VolumetricConstraints.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

std::shared_ptr<rigidbody::ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(const ConstraintManager* manager) {
    return generate_constraints(manager, settings::rigidbody::constraint_generation_strategy);
}

std::shared_ptr<rigidbody::ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(const ConstraintManager* manager, const settings::rigidbody::ConstraintGenerationStrategyChoice& choice) {
    switch (choice) {
        case settings::rigidbody::ConstraintGenerationStrategyChoice::Linear:
            return std::make_shared<rigidbody::LinearConstraints>(manager);
        case settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric:
            return std::make_shared<rigidbody::VolumetricConstraints>(manager);
        default: 
            throw except::unexpected("rigidbody::factory::generate_constraints: Unknown constraint generation strategy choice. Did you forget to add it to the switch statement?");
    }
}