#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <utility/Exceptions.h>

// extend the settings namespace with a new option
std::unique_ptr<rigidbody::TransformStrategy> rigidbody::factory::create_transform_strategy(rigidbody::RigidBody* body, settings::rigidbody::TransformationStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::TransformationStrategyChoice::RigidTransform:
            return std::make_unique<rigidbody::RigidTransform>(body); 
        case settings::rigidbody::TransformationStrategyChoice::SingleTransform:
            return std::make_unique<rigidbody::SingleTransform>(body);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_transform_strategy: Unkown TransformationStrategy. Did you forget to add it to the switch statement?");
    }
}