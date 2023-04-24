#include <settings/EMSettings.h>
#include <settings/SettingsIORegistry.h>

namespace settings::em {
    unsigned int sample_frequency = 1;
    double concentration = 1;
    unsigned int charge_levels = 1;
    bool hydrate = true;
    bool save_pdb = true;
    Limit alpha_levels = {0.5, 8};
    bool fixed_weights = false;
    bool plot_landscapes = false;
    bool simulation::noise = true;

    namespace io {
        settings::io::SettingSection general_settings("EM", {
            settings::io::create(sample_frequency, "sample_frequency"),
            settings::io::create(concentration, "concentration"),
            settings::io::create(charge_levels, "charge_levels"),
            settings::io::create(hydrate, "hydrate"),
            settings::io::create(save_pdb, "save_pdb"),
            settings::io::create(alpha_levels, "alpha_levels"),
            settings::io::create(fixed_weights, "fixed_weights"),
            settings::io::create(plot_landscapes, "plot_landscapes"),
            settings::io::create(simulation::noise, "simulation.noise")
        });
    }
}