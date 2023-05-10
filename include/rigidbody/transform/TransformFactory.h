#pragma once

#include <rigidbody/transform/TransformStrategy.h>

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a transformation strategy.
         */
        std::unique_ptr<TransformStrategy> create_transform_strategy(RigidBody* body);

        /**
         * @brief Prepare a transformation strategy.
         */
        std::unique_ptr<TransformStrategy> create_transform_strategy(RigidBody* body, const settings::rigidbody::TransformationStrategyChoice& choice);
    }
}