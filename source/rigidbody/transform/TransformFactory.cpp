#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <utility/Exceptions.h>

using namespace settings::rigidbody;
using namespace settings::detail;

// extend the settings namespace with a new option
template<> std::string SmartOption<TransformationStrategyChoice>::get() const {return std::to_string(static_cast<int>(value));}
template<> void SmartOption<TransformationStrategyChoice>::set(const std::vector<std::string>& val) {
    value = static_cast<TransformationStrategyChoice>(std::stoi(val[0]));
}

SmartOption<TransformationStrategyChoice> culling_strategy(TransformationStrategyChoice::RigidTransform);
std::unique_ptr<rigidbody::TransformStrategy> create_transform_strategy(RigidBody* body, TransformationStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::TransformationStrategyChoice::RigidTransform:
            return std::make_unique<rigidbody::RigidTransform>(body); 
        case settings::rigidbody::TransformationStrategyChoice::SingleTransform:
            return std::make_unique<rigidbody::SingleTransform>(body);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_transform_strategy: Unkown TransformationStrategy. Did you forget to add it to the switch statement?");
    }
}