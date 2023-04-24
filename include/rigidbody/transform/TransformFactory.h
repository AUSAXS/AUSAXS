#pragma once

#include <rigidbody/transform/TransformStrategy.h>
#include <settings/RigidBodySettings.h>

namespace rigidbody {
    namespace factory {
        /**
         * @brief Prepare a transformation strategy.
         */
        std::unique_ptr<TransformStrategy> create_transform_strategy(RigidBody* body, settings::rigidbody::TransformationStrategyChoice choice = settings::rigidbody::transform_strategy);
    }
}