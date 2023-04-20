#pragma once

#include <hydrate/placement/PlacementStrategy.h>

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

namespace grid {
    namespace factory {
        /**
         * @brief Prepare a placement class. 
         */
        std::unique_ptr<PlacementStrategy> construct_placement_strategy(Grid* grid, setting::grid::PlacementStrategy choice = setting::grid::placement_strategy);
    }
}