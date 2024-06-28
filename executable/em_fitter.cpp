#include <CLI/CLI.hpp>

#include <plots/All.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <utility/Console.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>
#include <constants/Constants.h>

#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    settings::em::mass_axis = true;
    settings::em::hydrate = true;
    settings::fit::verbose = true;
    settings::em::alpha_levels = {1, 10};
    settings::hist::weighted_bins = true;

    io::ExistingFile mfile, mapfile, settings;
    CLI::App app{"Fit an EM map to a SAXS measurement."};
    app.add_option("input-map", mapfile, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input-exp", mfile, "Path to the SAXS measurement.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/em_fitter/");
    app.add_option_function<std::string>("--unit,-u", [] (const std::string& s) {settings::detail::parse_option("unit", {s});}, "The unit of the q values in the measurement file. Options: A, nm.");
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads);
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
    app.add_option("--levelmin", settings::em::alpha_levels.min, "Lower limit on the alpha levels to use for the EM map. Note that lowering this limit severely impacts the performance and memory load.");
    app.add_option("--levelmax", settings::em::alpha_levels.max, "Upper limit on the alpha levels to use for the EM map. Increasing this limit improves the performance.");
    app.add_option("--charge-levels", settings::em::charge_levels, "Number of charge levels to use for the EM map.");
    app.add_option("--frequency", settings::em::sample_frequency, "Sampling frequency of the EM map.");
    app.add_option("--max-iterations", settings::fit::max_iterations, "Maximum number of iterations to perform. This is only approximate.");
    app.add_flag("--hydrate,!--no-hydrate", settings::em::hydrate, "Generate a hydration shell for the protein before fitting.");
    app.add_flag("--fixed-weight,!--dynamic-weight", settings::em::fixed_weights, "Use a fixed weight for the fit.");
    app.add_flag("--verbose,!--quiet", settings::fit::verbose, "Print the progress of the fit to the console.");
    app.add_flag("--weighted-bins, --!no-weighted-bins", settings::hist::weighted_bins, "Use weighted bins for the distance histograms.")->default_val(settings::hist::weighted_bins)->group("");
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");
    CLI11_PARSE(app, argc, argv);

    console::print_info("Running AUSAXS " + std::string(constants::version));

    //###################//
    //### PARSE INPUT ###//
    //###################//
    try {
        // if a settings file was provided
        if (p_settings->count() != 0) {
            settings::read(settings);        // read it
            CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
        } else {                            // otherwise check if there is a settings file in the same directory
            if (settings::discover(mfile.directory())) {
                CLI11_PARSE(app, argc, argv);
            }
        }
        if (!settings::general::output.empty() && settings::general::output.back() != '/') {settings::general::output += "/";}
        settings::general::output += mfile.stem() + "/" + mapfile.stem() + "/";

        // validate input
        if (!constants::filetypes::em_map.validate(mapfile)) {
            if (constants::filetypes::em_map.validate(mfile)) {
                std::swap(mapfile, mfile);
            } else {
                throw except::invalid_argument("Unknown EM extensions: \"" + mapfile + "\" and \"" + mfile + "\"");
            }
        }
        if (!constants::filetypes::saxs_data.validate(mfile)) {
            throw except::invalid_argument("Unknown SAXS data extension: \"" + mfile + "\"");
        }

        std::cout << "Performing EM fit with map \"" << mapfile << "\" and measurement \"" << mfile << "\"" << std::endl;
        em::ImageStack map(mapfile); 
        auto res = map.fit(mfile);

        fitter::FitReporter::report(res.get());
        fitter::FitReporter::save(res.get(), settings::general::output + "report.txt", argc, argv);

        res->info.dataset.save(settings::general::output + mfile.stem() + ".scat");
        res->info.fitted_intensity_interpolated.save(settings::general::output + "ausaxs.fit");
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        throw e;
    }
    return 0;
}