#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>
#include <utility/SmartOption.h>

namespace settings {
    namespace rigidbody {
        enum class ConstraintGenerationStrategyChoice {
            None,       // Do not generate constraints. Only those supplied by the user will be used.
            Linear,     // Generate a linear chain of constraints between bodies.
            Volumetric  // Generate constraints between bodies based on proximity. 
        };

        extern settings::detail::SmartOption<ConstraintGenerationStrategyChoice> constraint_generation_strategy;
    }
}

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::shared_ptr<ConstraintGenerationStrategy> generate_constraints(const ConstraintManager* manager, settings::rigidbody::ConstraintGenerationStrategyChoice choice = settings::rigidbody::constraint_generation_strategy);
    }
}