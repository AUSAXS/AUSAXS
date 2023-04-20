#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>
#include <utility/SmartOption.h>

#include <memory>

namespace settings {
    namespace rigidbody {
        enum class BodySelectStrategyChoice {
            RandomSelect,           // Select a random body, then a random constraint within that body. 
            RandomConstraintSelect, // Select a random constraint. 
            SequentialSelect        // Select the first constraint, then the second, etc.
        };

        extern settings::detail::SmartOption<BodySelectStrategyChoice> body_select_strategy;
    }
}

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a body selection strategy.
         */
        std::unique_ptr<BodySelectStrategy> create_selection_strategy(const RigidBody*, settings::rigidbody::BodySelectStrategyChoice choice = settings::rigidbody::body_select_strategy);
    }
}