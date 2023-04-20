#pragma once

#include <rigidbody/transform/TransformStrategy.h>
#include <utility/SmartOption.h>

namespace settings {
    namespace rigidbody {
        enum class TransformationStrategyChoice {
            RigidTransform,     // Transform all bodies connected to one side of the constraint. 
            SingleTransform,    // Transform only the body directly connected to one side of the constraint.
            ForceTransform      // Rotations and translations are applied as forces, resulting in more natural conformations. 
        };
        extern settings::detail::SmartOption<TransformationStrategyChoice> transform_strategy;
    }
}

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a transformation strategy.
         */
        std::unique_ptr<TransformStrategy> create_transform_strategy(RigidBody* body, settings::rigidbody::TransformationStrategyChoice choice = settings::rigidbody::transform_strategy);
    }
}