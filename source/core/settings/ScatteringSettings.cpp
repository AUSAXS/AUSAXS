// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/ScatteringSettings.h>

#include <form_factor/lookup/FormFactorManager.h>
#include <settings/SettingsIORegistry.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

settings::detail::Setting<settings::scattering::Radiation> settings::scattering::radiation{
    settings::scattering::Radiation::XRay,
    [] (settings::scattering::Radiation&) {
        ausaxs::form_factor::manager::rebuild(); // the form factor tables depend on the probe
    }
};

namespace {
    using namespace ausaxs::settings;
    settings::io::SettingSection scattering_section("Scattering", {
        settings::io::create(settings::scattering::radiation, "radiation")
    });
}

template<> std::string settings::io::detail::SettingRef<settings::scattering::Radiation>::get() const {
    switch (settingref) {
        case settings::scattering::Radiation::XRay:    return "xray";
        case settings::scattering::Radiation::Neutron: return "neutron";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::scattering::Radiation>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]);
    if (     str == "xray" || str == "x-ray") {settingref = settings::scattering::Radiation::XRay;}
    else if (str == "neutron") {settingref = settings::scattering::Radiation::Neutron;}
    else {
        throw except::io_error("settings: Unknown radiation \"" + str + "\". Valid options are \"xray\" and \"neutron\".");
    }
}
