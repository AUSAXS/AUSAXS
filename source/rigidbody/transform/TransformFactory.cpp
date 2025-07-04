// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;

std::unique_ptr<transform::TransformStrategy> rigidbody::factory::create_transform_strategy(observer_ptr<rigidbody::RigidBody> body) {
    return create_transform_strategy(body, settings::rigidbody::transform_strategy);
}

std::unique_ptr<transform::TransformStrategy> rigidbody::factory::create_transform_strategy(
    observer_ptr<rigidbody::RigidBody> body, const settings::rigidbody::TransformationStrategyChoice& choice
) {
    switch (choice) {
        case settings::rigidbody::TransformationStrategyChoice::RigidTransform:
            return std::make_unique<transform::RigidTransform>(body); 
        case settings::rigidbody::TransformationStrategyChoice::SingleTransform:
            return std::make_unique<transform::SingleTransform>(body);
        default: 
            throw except::unknown_argument(
                "rigidbody::factory::create_transform_strategy: Unkown TransformationStrategy. "
                "Did you forget to add it to the switch statement?"
            );
    }
}