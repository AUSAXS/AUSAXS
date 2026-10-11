// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/EMSettings.h>

#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

settings::Setting<int> settings::em::sample_frequency{1};
settings::Setting<double> settings::em::concentration{1};
settings::Setting<int> settings::em::charge_levels{50};
settings::Setting<bool> settings::em::hydrate{true};
settings::Setting<bool> settings::em::save_pdb{true};
settings::Setting<Limit> settings::em::alpha_levels{{1, 10}};
settings::Setting<bool> settings::em::fixed_weights{true};
settings::Setting<bool> settings::em::plot_landscapes{false};
settings::Setting<bool> settings::em::simulation::noise{true};
settings::Setting<bool> settings::em::mass_axis{true};

namespace {
    using namespace ausaxs::settings;
    settings::io::SettingSection em_section("EM", {
        settings::io::create(em::sample_frequency, "sample_frequency"),
        settings::io::create(em::concentration, "concentration"),
        settings::io::create(em::charge_levels, "charge_levels"),
        settings::io::create(em::hydrate, "hydrate"),
        settings::io::create(em::mass_axis, "mass_axis"),
        settings::io::create(em::save_pdb, "save_pdb"),
        settings::io::create(em::alpha_levels, "alpha_levels"),
        settings::io::create(em::fixed_weights, "fixed_weights"),
        settings::io::create(em::plot_landscapes, "plot_landscapes"),
        settings::io::create(em::simulation::noise, "simulation.noise")
    });
}