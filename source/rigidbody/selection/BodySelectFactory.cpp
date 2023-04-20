#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/selection/RandomSelect.h>
#include <rigidbody/selection/SequentialSelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <utility/Exceptions.h>

using namespace settings::rigidbody;
using namespace settings::detail;

// extend the settings namespace with a new option
template<> std::string SmartOption<BodySelectStrategyChoice>::get() const {return std::to_string(static_cast<int>(value));}
template<> void SmartOption<BodySelectStrategyChoice>::set(const std::vector<std::string>& val) {
    value = static_cast<BodySelectStrategyChoice>(std::stoi(val[0]));
}

SmartOption<BodySelectStrategyChoice> culling_strategy(BodySelectStrategyChoice::RandomSelect);
std::unique_ptr<rigidbody::BodySelectStrategy> create_selection_strategy(const RigidBody* body, BodySelectStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::BodySelectStrategyChoice::RandomSelect:
            return std::make_unique<rigidbody::RandomSelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect:
            return std::make_unique<rigidbody::RandomConstraintSelect>(body);
        case settings::rigidbody::BodySelectStrategyChoice::SequentialSelect:
            return std::make_unique<rigidbody::SequentialSelect>(body);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_selection_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}