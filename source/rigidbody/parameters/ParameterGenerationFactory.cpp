#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/SimpleParameterGeneration.h>
#include <rigidbody/parameters/RotationsOnly.h>
#include <utility/Exceptions.h>

using namespace settings::rigidbody;
using namespace settings::detail;

// extend the settings namespace with a new option
template<> std::string SmartOption<ParameterGenerationStrategyChoice>::get() const {return std::to_string(static_cast<int>(value));}
template<> void SmartOption<ParameterGenerationStrategyChoice>::set(const std::vector<std::string>& val) {
    value = static_cast<ParameterGenerationStrategyChoice>(std::stoi(val[0]));
}

SmartOption<ParameterGenerationStrategyChoice> culling_strategy(ParameterGenerationStrategyChoice::Simple);
std::unique_ptr<rigidbody::ParameterGenerationStrategy> create_parameter_strategy(unsigned int iterations, double translate_amp, double rotate_amp, ParameterGenerationStrategyChoice choice) {
    switch (choice) {
        case settings::rigidbody::ParameterGenerationStrategyChoice::Simple:
            return std::make_unique<rigidbody::SimpleParameterGeneration>(iterations, translate_amp, rotate_amp);
        case settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly:
            return std::make_unique<rigidbody::RotationsOnly>(iterations, translate_amp, rotate_amp);
        default: 
            throw except::unknown_argument("rigidbody::factory::create_parameter_strategy: Unknown strategy. Did you forget to add it to the switch statement?");
    }
}