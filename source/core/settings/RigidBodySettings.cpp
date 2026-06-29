// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/RigidBodySettings.h>
#include <settings/SettingRef.h>
#include <settings/SettingsIORegistry.h>
#include <settings/MoleculeSettings.h>

using namespace ausaxs;

unsigned int settings::rigidbody::iterations = 1000;
double settings::rigidbody::bond_distance = 3;
settings::rigidbody::TransformationStrategyChoice settings::rigidbody::transform_strategy = TransformationStrategyChoice::RigidTransform;
settings::rigidbody::ParameterGenerationStrategyChoice settings::rigidbody::parameter_generation_strategy = ParameterGenerationStrategyChoice::Simple;
settings::rigidbody::ParameterMaskStrategyChoice settings::rigidbody::parameter_mask_strategy = ParameterMaskStrategyChoice::All;
settings::rigidbody::BodySelectStrategyChoice settings::rigidbody::body_select_strategy = BodySelectStrategyChoice::RandomBodySelect;
settings::detail::Setting<settings::rigidbody::ConstraintGenerationStrategyChoice> settings::rigidbody::constraint_generation_strategy = {
    settings::rigidbody::ConstraintGenerationStrategyChoice::None,
    [] (settings::rigidbody::ConstraintGenerationStrategyChoice& val) {
        // backbone-based constraint generators rely on per-atom C-alpha metadata, which is only
        // retained at load time when store_calpha is enabled; couple the two so callers cannot forget.
        // residue_seq is likewise enabled so the generators can prefer sequential C-alpha pairs when
        // possible; it remains a soft preference (the generators degrade gracefully if it is absent).
        if (val == settings::rigidbody::ConstraintGenerationStrategyChoice::Linear ||
            val == settings::rigidbody::ConstraintGenerationStrategyChoice::Volumetric)
        {
            settings::molecule::store_calpha = true;
            settings::molecule::store_residue_seq = true;
        }
    }
};
settings::rigidbody::DecayStrategyChoice settings::rigidbody::decay_strategy = DecayStrategyChoice::Linear;
settings::rigidbody::ControllerChoice settings::rigidbody::controller_choice = ControllerChoice::Classic;

std::vector<int> ausaxs::settings::rigidbody::detail::constraints;
std::string ausaxs::settings::rigidbody::detail::calibration_file;

namespace ausaxs::settings::io {
    settings::io::SettingSection rigidbody_section("RigidBody", {
        settings::io::create(rigidbody::iterations, "iterations"),
        settings::io::create(rigidbody::bond_distance, "bond_distance"),
        settings::io::create(rigidbody::detail::constraints, "constraints"),
        settings::io::create(rigidbody::detail::calibration_file, "calibration_file")
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

template<> std::string settings::io::detail::SettingRef<settings::rigidbody::DecayStrategyChoice>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::rigidbody::DecayStrategyChoice>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::rigidbody::DecayStrategyChoice>(std::stoi(val[0]));
}

template<> std::string settings::io::detail::SettingRef<settings::rigidbody::ControllerChoice>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::rigidbody::ControllerChoice>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::rigidbody::ControllerChoice>(std::stoi(val[0]));
}