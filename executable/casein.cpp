#include <CLI/CLI.hpp>

#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <data/Protein.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>
#include <utility/Constants.h>

int main(int argc, char const *argv[]) {
    settings::protein::use_effective_charge = false;
    settings::fit::verbose = true;

    CLI::App app{"Casein micelle analysis."};
    io::ExistingFile mapfile, settings;
    app.add_option("input-map", mapfile, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/casein/");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
    app.add_option("--frequency", settings::em::sample_frequency, "Sampling frequency of the EM map.");
    app.add_flag("--hydrate,!--no-hydrate", settings::em::hydrate, "Whether to hydrate the protein before fitting.");
    app.add_flag("--fixed-weight,!--no-fixed-weight", settings::em::fixed_weights, "Whether to use a fixed weight for the fit.");
    CLI11_PARSE(app, argc, argv);

    //###################//
    //### PARSE INPUT ###//
    //###################//
    settings::general::output += mapfile.stem() + "/";

    // if a settings file was provided
    if (p_settings->count() != 0) {
        settings::read(settings);       // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (settings::discover(std::filesystem::path(mapfile).parent_path().string())) {
            CLI11_PARSE(app, argc, argv);
        }
    }

    em::ImageStack map(mapfile);
    auto protein = map.get_protein(10);
    plots::PlotIntensity(protein->get_histogram(), settings::general::output + "intensity.png");
}