// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/ScatteringSettings.h>

#include <constants/Constants.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/neutron/NeutronSolvent.h>
#include <settings/SettingsIORegistry.h>
#include <utility/Console.h>
#include <utility/StringUtils.h>

#include <string>

using namespace ausaxs;

settings::detail::Setting<settings::scattering::Radiation> settings::scattering::radiation{
    settings::scattering::Radiation::XRay,
    [] (settings::scattering::Radiation&) {
        ausaxs::form_factor::manager::rebuild(); // the form factor tables depend on the probe
    }
};

settings::detail::Setting<double> settings::scattering::xray::solvent_density{
    constants::charge::density::water,
    [] (double&) {
        ausaxs::form_factor::manager::rebuild(); // the excluded volume tables are scaled by the solvent density
    }
};

settings::detail::Setting<double> settings::scattering::neutron::solvent_density{
    ausaxs::form_factor::neutron::solvent_density(1),
    [] (double& density) {
        if (double x = ausaxs::form_factor::neutron::d2o_fraction(density); x < 0 || 1 < x) {
            console::print_warning(
                "settings::scattering::neutron::solvent_density: The density " + std::to_string(density) + " fm/Å^3 lies outside the range of H2O/D2O mixtures "
                "[" + std::to_string(ausaxs::form_factor::neutron::solvent_density(0)) + ", " + std::to_string(ausaxs::form_factor::neutron::solvent_density(1)) + "]. "
                "The isotope composition of the hydration water and the exchanged hydrogens is extrapolated to a D2O fraction of " + std::to_string(x) + "."
            );
        }
        ausaxs::form_factor::manager::rebuild(); // the excluded volume tables are scaled by the solvent density, and the hydrogens take its isotope composition
    }
};

settings::detail::Setting<double> settings::scattering::neutron::exchange_efficiency{
    1,
    [] (double&) {
        ausaxs::form_factor::manager::rebuild(); // the exchangeable hydrogens take the isotope composition of the solvent
    }
};

void settings::scattering::neutron::set_d2o_fraction(double d2o_fraction) {
    solvent_density = ausaxs::form_factor::neutron::solvent_density(d2o_fraction);
}

namespace {
    using namespace ausaxs::settings;
    settings::io::SettingSection scattering_section("Scattering", {
        settings::io::create(settings::scattering::radiation, "radiation"),
        settings::io::create(settings::scattering::xray::solvent_density, "xray_solvent_density"),
        settings::io::create(settings::scattering::neutron::solvent_density, "neutron_solvent_density"),
        settings::io::create(settings::scattering::neutron::exchange_efficiency, "neutron_exchange_efficiency")
    });
}

template<> std::string settings::io::detail::SettingRef<settings::scattering::Radiation>::get() const {
    switch (settingref) {
        case settings::scattering::Radiation::XRay:    return "xray";
        case settings::scattering::Radiation::Neutron: return "neutron";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> std::string settings::io::detail::SettingRef<settings::scattering::Radiation>::type() const {return "string";}

template<> void settings::io::detail::SettingRef<settings::scattering::Radiation>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]);
    if (     str == "xray" || str == "x-ray") {settingref = settings::scattering::Radiation::XRay;}
    else if (str == "neutron") {settingref = settings::scattering::Radiation::Neutron;}
    else {
        throw except::io_error("settings: Unknown radiation \"" + str + R"(". Valid options are "xray" and "neutron".)");
    }
}
