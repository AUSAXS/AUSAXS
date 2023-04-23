#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>
#include <rigidbody/RigidBodySettings.h>

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::shared_ptr<ConstraintGenerationStrategy> generate_constraints(const ConstraintManager* manager, settings::rigidbody::ConstraintGenerationStrategyChoice choice = settings::rigidbody::constraint_generation_strategy);
    }
}