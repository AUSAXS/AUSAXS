#include <rigidbody/constraints/generation/ConstraintGenerationFactory.h>
#include <rigidbody/constraints/generation/LinearConstraints.h>
#include <rigidbody/constraints/generation/VolumetricConstraints.h>
#include <utility/Exceptions.h>

using namespace settings::rigidbody;
using namespace settings::detail;

// extend the settings namespace with a new option
template<> std::string SmartOption<ConstraintGenerationStrategyChoice>::get() const {return std::to_string(static_cast<int>(value));}
template<> void SmartOption<ConstraintGenerationStrategyChoice>::set(const std::vector<std::string>& val) {
    value = static_cast<ConstraintGenerationStrategyChoice>(std::stoi(val[0]));
}

SmartOption<ConstraintGenerationStrategyChoice> culling_strategy(ConstraintGenerationStrategyChoice::Linear);
std::shared_ptr<rigidbody::ConstraintGenerationStrategy> rigidbody::factory::generate_constraints(const ConstraintManager* manager, ConstraintGenerationStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::ConstraintGenerationStrategyChoice::Linear:
            return std::make_shared<rigidbody::LinearConstraints>(manager);
        case settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric:
            return std::make_shared<rigidbody::VolumetricConstraints>(manager);
        default: 
            throw except::unexpected("rigidbody::factory::generate_constraints: Unknown constraint generation strategy choice. Did you forget to add it to the switch statement?");
    }
}