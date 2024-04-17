/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/RigidBodySettings.h>
#include <settings/SettingRef.h>
#include <settings/SettingsIORegistry.h>

unsigned int settings::rigidbody::iterations = 1000;
double settings::rigidbody::bond_distance = 3;
settings::rigidbody::TransformationStrategyChoice settings::rigidbody::transform_strategy = TransformationStrategyChoice::RigidTransform;
settings::rigidbody::ParameterGenerationStrategyChoice settings::rigidbody::parameter_generation_strategy = ParameterGenerationStrategyChoice::Simple;
settings::rigidbody::BodySelectStrategyChoice settings::rigidbody::body_select_strategy = BodySelectStrategyChoice::RandomBodySelect;
settings::rigidbody::ConstraintGenerationStrategyChoice settings::rigidbody::constraint_generation_strategy = ConstraintGenerationStrategyChoice::Linear;
settings::rigidbody::DecayStrategyChoice settings::rigidbody::decay_strategy = DecayStrategyChoice::Linear;

namespace settings::rigidbody::detail {
    std::vector<int> constraints;
    std::string calibration_file;
}

namespace settings::rigidbody::io {
    settings::io::SettingSection rigidbody_settings("RigidBody", {
        settings::io::create(iterations, "iterations"),
        settings::io::create(bond_distance, "bond_distance"),
        settings::io::create(detail::constraints, "constraints"),
        settings::io::create(detail::calibration_file, "calibration_file")
    });
}

template<> std::string settings::io::detail::SettingRef<settings::rigidbody::TransformationStrategyChoice>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::rigidbody::TransformationStrategyChoice>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::rigidbody::TransformationStrategyChoice>(std::stoi(val[0]));
}

template<> std::string settings::io::detail::SettingRef<settings::rigidbody::ParameterGenerationStrategyChoice>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::rigidbody::ParameterGenerationStrategyChoice>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::rigidbody::ParameterGenerationStrategyChoice>(std::stoi(val[0]));
}

template<> std::string settings::io::detail::SettingRef<settings::rigidbody::BodySelectStrategyChoice>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::rigidbody::BodySelectStrategyChoice>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::rigidbody::BodySelectStrategyChoice>(std::stoi(val[0]));
}

template<> std::string settings::io::detail::SettingRef<settings::rigidbody::ConstraintGenerationStrategyChoice>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::rigidbody::ConstraintGenerationStrategyChoice>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::rigidbody::ConstraintGenerationStrategyChoice>(std::stoi(val[0]));
}