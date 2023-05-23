#pragma once

#include <rigidbody/constraints/generation/ConstraintGenerationStrategy.h>

namespace settings::rigidbody {enum class ConstraintGenerationStrategyChoice;}
namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<ConstraintGenerationStrategy> generate_constraints(const ConstraintManager* manager);

        /**
         * @brief Prepare a constraint generator. 
         */
        std::unique_ptr<ConstraintGenerationStrategy> generate_constraints(const ConstraintManager* manager, const settings::rigidbody::ConstraintGenerationStrategyChoice& choice);
    }
}