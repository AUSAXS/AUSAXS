#include <CLI/CLI.hpp>

#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::hydrate = true;
    setting::fit::verbose = true;
    // setting::em::alpha_levels = {0.05, 8};

    CLI::App app{"Fit an EM map to a SAXS measurement."};

    std::string mfile, mapfile, settings, output;
    app.add_option("input-map", mapfile, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input-exp", mfile, "Path to the SAXS measurement.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", output, "Path to save the generated figures at.");
    app.add_option("--qmin", setting::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", setting::axes::qmax, "Upper limit on used q values from measurement file.");
    app.add_option("--levelmin", setting::em::alpha_levels.min, "Lower limit on the alpha levels to use for the EM map. Note that lowering this limit severely impacts the performance.");
    app.add_option("--levelmax", setting::em::alpha_levels.max, "Upper limit on the alpha levels to use for the EM map.");
    app.add_option("--frequency", setting::em::sample_frequency, "Sampling frequency of the EM map.");
    app.add_flag("--hydrate,!--no-hydrate", setting::em::hydrate, "Whether to hydrate the protein before fitting.");
    app.add_flag("--fixed-weight,!--no-fixed-weight", setting::em::fixed_weights, "Whether to use a fixed weight for the fit.");
    CLI11_PARSE(app, argc, argv);

    // if a settings file was provided
    if (p_settings->count() != 0) {
        setting::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        setting::discover(std::filesystem::path(mfile).parent_path().string());
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

    if (output.empty()) {
        output = "figures/";
    }
    setting::plot::path = output + "em_fitter/" + utility::stem(mfile) + "/" + utility::stem(mapfile) + "/";
    
    std::cout << "Performing EM fit with map " << mapfile << " and measurement " << mfile << std::endl;
    em::ImageStack map(mapfile); 

    // Fit the measurements to the EM density map.
    auto res = map.fit(mfile);
    FitReporter::report(res);
    FitReporter::save(res, setting::plot::path + "report.txt");

    res->figures.data.save(setting::plot::path + utility::stem(mfile) + ".dat");
    res->figures.intensity_interpolated.save(setting::plot::path + "fit.fit");
    plots::PlotIntensityFit::quick_plot(res, setting::plot::path + "intensity_fit." + setting::plot::format);
    plots::PlotIntensityFitResiduals::quick_plot(res, setting::plot::path + "residuals." + setting::plot::format);
    return 0;
}