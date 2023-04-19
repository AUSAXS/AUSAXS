#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <utility/Exceptions.h>

std::shared_ptr<rigidbody::ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(const rigidbody::ConstraintManager* manager, setting::rigidbody::ConstraintGenerationStrategyChoice choice) {
    switch (choice) {
        case setting::rigidbody::ConstraintGenerationStrategyChoice::Linear:
            return std::make_shared<rigidbody::LinearConstraints>(manager);
        case setting::rigidbody::ConstraintGenerationStrategyChoice::Volumetric:
            return std::make_shared<rigidbody::VolumetricConstraints>(manager);
        default: 
            throw except::unexpected("rigidbody::factory::generate_constraints: Unknown constraint generation strategy choice. Did you forget to add it to the switch statement?");
    }
}