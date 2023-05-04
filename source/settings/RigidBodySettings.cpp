#include <settings/RigidBodySettings.h>
#include <settings/SettingRef.h>
#include <settings/SettingsIORegistry.h>

namespace settings::rigidbody {
    unsigned int iterations = 1000;
    double bond_distance = 3;
    TransformationStrategyChoice transform_strategy = TransformationStrategyChoice::RigidTransform;
    ParameterGenerationStrategyChoice parameter_generation_strategy = ParameterGenerationStrategyChoice::Simple;
    BodySelectStrategyChoice body_select_strategy = BodySelectStrategyChoice::RandomSelect;
    ConstraintGenerationStrategyChoice constraint_generation_strategy = ConstraintGenerationStrategyChoice::Linear;
    DecayStrategyChoice decay_strategy = DecayStrategyChoice::Linear;

    namespace detail {
        std::vector<int> constraints;
        std::string calibration_file;
    }

    namespace io {
        settings::io::SettingSection rigidbody_settings("RigidBody", {
            settings::io::create(iterations, "iterations"),
            settings::io::create(bond_distance, "bond_distance"),
            settings::io::create(detail::constraints, "constraints"),
            settings::io::create(detail::calibration_file, "calibration_file")
        });
    }
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