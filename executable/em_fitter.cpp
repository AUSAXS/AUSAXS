#include <CLI/CLI.hpp>

#include <plots/All.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>
#include <constants/Constants.h>

#include <iostream>

using namespace ausaxs;

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    settings::fit::verbose = true;
    settings::general::supplementary_plots = false;
    bool save_settings = false;

    io::ExistingFile mfile, mapfile, settings;
    CLI::App app{"Fit an EM map to a SAXS measurement."};
    auto input_map = app.add_option("input-map", mapfile, "Path to the EM map.")->check(CLI::ExistingFile);
    auto input_saxs = app.add_option("input-saxs", mfile, "Path to the SAXS measurement.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Output folder to write the results to.")->default_val("output/em_fitter/");
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");
    app.add_flag_callback("-v,--version", [] () {std::cout << constants::version << std::endl; exit(0);}, "Print the AUSAXS version.");
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads);

    // config subcommands
    auto sub_config = app.add_subcommand("config", "See and set additional options for the configuration.");
    auto p_settings = sub_config->add_option("--file,-f", settings, "The configuration file to use.")->check(CLI::ExistingFile);
    sub_config->add_flag("--save", save_settings, "Save the settings to a file.");
    sub_config->add_flag_callback("--log", [] () {logging::start("em_fitter");}, "Enable logging to a file.");

    // data subcommands
    auto sub_data = app.add_subcommand("saxs", "See and set additional options for the SAXS data.");
    sub_data->add_option(
        "--qmax", 
        settings::axes::qmax, 
        "Upper limit on used q values from the measurement file.")
        ->default_val(settings::axes::qmax)
        ->check(CLI::Range(constants::axes::q_axis.min, constants::axes::q_axis.max))
    ;
    sub_data->add_option(
        "--qmin", 
        settings::axes::qmin, 
        "Lower limit on used q values from the measurement file.")
        ->default_val(settings::axes::qmin)
        ->check(CLI::Range(constants::axes::q_axis.min, constants::axes::q_axis.max))
    ;
    sub_data->add_option_function<std::string>("--unit,-u", [] (const std::string& s) {settings::detail::parse_option("unit", {s});}, "The unit of the q values in the measurement file. Options: A, nm.");
    sub_data->add_option("--skip", settings::axes::skip, "Number of points to skip in the measurement file.")->default_val(settings::axes::skip);
    sub_data->add_flag("--rebin", settings::flags::data_rebin, "Rebin the data to increase the information content of each data point.")->default_val(settings::flags::data_rebin);

    // em subcommands
    auto sub_em = app.add_subcommand("em", "See and set additional options for the EM map.");
    sub_em->add_option("--levelmin", settings::em::alpha_levels.min, "Lower limit on the alpha levels to use for the EM map. Note that lowering this limit severely impacts the performance and memory load.");
    sub_em->add_option("--levelmax", settings::em::alpha_levels.max, "Upper limit on the alpha levels to use for the EM map. Increasing this limit improves the performance.");
    sub_em->add_option("--charge-levels", settings::em::charge_levels, "Number of charge levels to use for the EM map.");
    sub_em->add_option("--frequency", settings::em::sample_frequency, "Sampling frequency of the EM map.");
    sub_em->add_flag("--hydrate,!--no-hydrate", settings::em::hydrate, "Generate a hydration shell for the protein before fitting.");
    sub_em->add_flag("--fixed-weight,!--dynamic-weight", settings::em::fixed_weights, "Use a fixed weight for the fit.");

    // fit subcommands
    auto sub_fit = app.add_subcommand("fit", "See and set additional options for the fitting process.");
    sub_fit->add_option("--max-iterations", settings::fit::max_iterations, "Maximum number of iterations to perform. This is only approximate.");
    sub_fit->add_flag("--verbose,!--quiet", settings::fit::verbose, "Print the progress of the fit to the console.");

    app.final_callback([&] () {
        // save settings if requested
        if (save_settings) {
            settings::write("settings.txt");
            console::print_info("Settings saved to settings.txt in current directory.");
            if (!input_map->count() || !input_saxs->count()) { // gracefully exit if no input files are provided
                exit(0);
            }
        }

        // required args (not marked ->required() since that interferes with the help flag for subcommands)
        if (!input_map->count() || !input_saxs->count()) {
            std::cout << "Error: Both input_structure and input_measurement are required." << std::endl;
            exit(1);
        }
    });

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

        // validate input
        if (!constants::filetypes::em_map.check(mapfile)) {
            if (constants::filetypes::em_map.check(mfile)) {
                std::swap(mapfile, mfile);
            } else {
                throw except::invalid_argument("Unknown EM extensions: \"" + mapfile.str() + "\" and \"" + mfile.str() + "\"");
            }
        }
        if (!constants::filetypes::saxs_data.check(mfile)) {
            throw except::invalid_argument("Unknown SAXS data extension: \"" + mfile.str() + "\"");
        }
        if (!settings::em::hydrate) {settings::fit::fit_hydration = false;} 
        if (!settings::general::output.empty() && settings::general::output.back() != '/') {settings::general::output += "/";}
        settings::general::output += mfile.stem() + "/" + mapfile.stem() + "/";

        std::cout << "Performing EM fit with map \"" << mapfile << "\" and measurement \"" << mfile << "\"" << std::endl;
        em::ImageStack map(mapfile); 
        auto res = map.fit(mfile);

        fitter::FitReporter::report(res.get());
        fitter::FitReporter::save(res.get(), {settings::general::output + "report.txt"}, argc, argv);
        res->curves.select_columns({0, 1, 2, 3}).save(
            settings::general::output + "ausaxs.fit",
            "chi2=" + std::to_string(res->fval/res->dof) + " dof=" + std::to_string(res->dof)
        );
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        throw e;
    }
    return 0;
}