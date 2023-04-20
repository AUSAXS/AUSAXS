#include <em/EMSettings.h>

namespace settings::em {
    settings::detail::SmartOption<unsigned int> sample_frequency(1, "sample-frequency");
    settings::detail::SmartOption<double> concentration(1, "concentration");
    settings::detail::SmartOption<unsigned int> charge_levels(1, "charge-levels");

    settings::detail::SmartOption<bool> hydrate(true, "hydrate");

    settings::detail::SmartOption<bool> save_pdb(true, "save-pdb");
    settings::detail::SmartOption<Limit> alpha_levels({0.5, 8}, {"alpha-levels", "alpha"});

    settings::detail::SmartOption<bool> fixed_weights(false, {"fixed-weights", "fixed-weight"});
    settings::detail::SmartOption<bool> plot_landscapes(false, "plot-landscapes");

    settings::detail::SmartOption<bool> simulation::noise(true, "noise");
}