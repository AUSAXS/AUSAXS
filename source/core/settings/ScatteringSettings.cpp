// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/ScatteringSettings.h>

#include <constants/Constants.h>
#include <form_factor/NeutronFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <settings/InternalState.h>
#include <settings/SettingsIORegistry.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

namespace {
    // Recompute the cached solvent density of the active probe. Must run before the form factor tables are rebuilt.
    void update_solvent_density() {
        using Radiation = settings::scattering::Radiation;
        switch (settings::scattering::radiation.value) {
            case Radiation::XRay:    settings::internal_state::solvent_density = settings::scattering::xray_solvent_density.value; break;
            case Radiation::Neutron: settings::internal_state::solvent_density = settings::scattering::neutron_solvent_density.value; break;
            default: throw except::unexpected("settings::scattering::update_solvent_density: Unknown radiation type.");
        }
        ausaxs::form_factor::ExvTableManager::clear_exv_form_factor_sets(); // the cached sets are built with the old density
    }
}

settings::detail::Setting<settings::scattering::Radiation> settings::scattering::radiation{
    settings::scattering::Radiation::XRay,
    [] (settings::scattering::Radiation&) {
        update_solvent_density();
        ausaxs::form_factor::manager::rebuild(); // the form factor tables depend on the probe
    }
};

settings::detail::Setting<double> settings::scattering::xray_solvent_density{
    constants::charge::density::water,
    [] (double&) {
        update_solvent_density();
        ausaxs::form_factor::manager::rebuild(); // the excluded volume tables are scaled by the solvent density
    }
};

settings::detail::Setting<double> settings::scattering::neutron_solvent_density{
    form_factor::neutron::solvent_density(0),
    [] (double&) {
        update_solvent_density();
        ausaxs::form_factor::manager::rebuild(); // the excluded volume tables are scaled by the solvent density
    }
};

namespace {
    using namespace ausaxs::settings;
    settings::io::SettingSection scattering_section("Scattering", {
        settings::io::create(settings::scattering::radiation, "radiation"),
        settings::io::create(settings::scattering::xray_solvent_density, "xray_solvent_density"),
        settings::io::create(settings::scattering::neutron_solvent_density, "neutron_solvent_density")
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
