/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <settings/RigidBodySettings.h>
#include <utility/Exceptions.h>

using namespace rigidbody;

std::unique_ptr<transform::TransformStrategy> rigidbody::factory::create_transform_strategy(rigidbody::RigidBody* body) {
    return create_transform_strategy(body, settings::rigidbody::transform_strategy);
}

std::unique_ptr<transform::TransformStrategy> rigidbody::factory::create_transform_strategy(rigidbody::RigidBody* body, const settings::rigidbody::TransformationStrategyChoice& choice) {
    switch (choice) {
        case settings::rigidbody::TransformationStrategyChoice::RigidTransform:
            return std::make_unique<transform::RigidTransform>(body); 
        case settings::rigidbody::TransformationStrategyChoice::SingleTransform:
            return std::make_unique<transform::SingleTransform>(body);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_transform_strategy: Unkown TransformationStrategy. Did you forget to add it to the switch statement?");
    }
}