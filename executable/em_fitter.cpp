#include <CLI/CLI.hpp>

#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>
#include <utility/Constants.h>

int main(int argc, char const *argv[]) {
    settings::protein::use_effective_charge = false;
    settings::em::mass_axis = true;
    settings::em::hydrate = true;
    settings::fit::verbose = true;
    settings::em::alpha_levels = {0.05, 5};

    io::ExistingFile mfile, mapfile, settings;
    CLI::App app{"Fit an EM map to a SAXS measurement."};
    app.add_option("input-map", mapfile, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input-exp", mfile, "Path to the SAXS measurement.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/em_fitter/");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
    app.add_option("--levelmin", settings::em::alpha_levels.min, "Lower limit on the alpha levels to use for the EM map. Note that lowering this limit severely impacts the performance.");
    app.add_option("--levelmax", settings::em::alpha_levels.max, "Upper limit on the alpha levels to use for the EM map. Increasing this limit improves the performance.");
    app.add_option("--frequency", settings::em::sample_frequency, "Sampling frequency of the EM map.");
    app.add_flag("--mass-axis,!--no-mass-axis", settings::em::mass_axis, "Whether to use a mass axis in place of the threshold axis.");
    app.add_flag("--hydrate,!--no-hydrate", settings::em::hydrate, "Whether to hydrate the protein before fitting.");
    app.add_flag("--fixed-weight,!--no-fixed-weight", settings::em::fixed_weights, "Whether to use a fixed weight for the fit.");
    CLI11_PARSE(app, argc, argv);

    //###################//
    //### PARSE INPUT ###//
    //###################//
    settings::general::output += mfile.stem() + "/" + mapfile.stem() + "/";

    // if a settings file was provided
    if (p_settings->count() != 0) {
        settings::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (settings::discover(std::filesystem::path(mfile).parent_path().string())) {
            CLI11_PARSE(app, argc, argv);
        }
    }

    // validate input
    if (!constants::filetypes::em_map.validate(mapfile)) {
        if (constants::filetypes::em_map.validate(mfile)) {
            std::swap(mapfile, mfile);
        } else {
            throw except::invalid_argument("Unknown EM extensions: " + mapfile + " and " + mfile);
        }
    }
    if (!constants::filetypes::saxs_data.validate(mfile)) {
        throw except::invalid_argument("Unknown SAXS data extension: " + mfile);
    }

    std::cout << "Performing EM fit with map " << mapfile << " and measurement " << mfile << std::endl;
    em::ImageStack map(mapfile); 

    // Fit the measurements to the EM density map.
    auto res = map.fit(mfile);

    std::string cmd_line;
    for (int i = 0; i < argc; ++i) {cmd_line.append(argv[i]).append(" ");}

    fitter::FitReporter::report(res);
    fitter::FitReporter::save(res, settings::general::output + "report.txt", cmd_line);

    res->figures.data.save(settings::general::output + mfile.stem() + ".dat");
    res->figures.intensity_interpolated.save(settings::general::output + "fit.fit");
    plots::PlotIntensityFit::quick_plot(res, settings::general::output + "intensity_fit." + settings::plots::format);
    plots::PlotIntensityFitResiduals::quick_plot(res, settings::general::output + "residuals." + settings::plots::format);
    return 0;
}