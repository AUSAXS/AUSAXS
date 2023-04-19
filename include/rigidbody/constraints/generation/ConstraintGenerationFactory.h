#pragma once

#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/generation/VolumetricConstraints.h>

//! REFACTOR
#include <utility/Settings.h>

// namespace setting {
//     namespace rigidbody {
//         enum class ConstraintGenerationStrategyChoice {
//             None,       // Do not generate constraints. Only those supplied by the user will be used.
//             Linear,     // Generate a linear chain of constraints between bodies.
//             Volumetric  // Generate constraints between bodies based on proximity. 
//         };

//         inline static ConstraintGenerationStrategyChoice csc = ConstraintGenerationStrategyChoice::Linear;
//     }
// }

namespace rigidbody {
    class ConstraintManager;

    namespace factory {
        /**
         * @brief Prepare a constraint generator. 
         */
        std::shared_ptr<ConstraintGenerationStrategy> generate_constraints(const ConstraintManager* manager, setting::rigidbody::ConstraintGenerationStrategyChoice choice = setting::rigidbody::csc);
    }
}