/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/EMSettings.h>
#include <settings/SettingsIORegistry.h>

unsigned int settings::em::sample_frequency = 1;
double settings::em::concentration = 1;
unsigned int settings::em::charge_levels = 50;
bool settings::em::hydrate = true;
bool settings::em::save_pdb = true;
Limit settings::em::alpha_levels = {0.5, 8};
bool settings::em::fixed_weights = true;
bool settings::em::plot_landscapes = false;
bool settings::em::simulation::noise = true;
bool settings::em::mass_axis = true;

namespace settings::em::io {
    settings::io::SettingSection general_settings("EM", {
        settings::io::create(sample_frequency, "sample_frequency"),
        settings::io::create(concentration, "concentration"),
        settings::io::create(charge_levels, "charge_levels"),
        settings::io::create(hydrate, "hydrate"),
        settings::io::create(mass_axis, "mass_axis"),
        settings::io::create(save_pdb, "save_pdb"),
        settings::io::create(alpha_levels, "alpha_levels"),
        settings::io::create(fixed_weights, "fixed_weights"),
        settings::io::create(plot_landscapes, "plot_landscapes"),
        settings::io::create(simulation::noise, "simulation.noise")
    });
}