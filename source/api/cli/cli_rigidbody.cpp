// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_rigidbody.h>
#include <CLI/CLI.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <fitter/FitReporter.h>
#include <constants/Constants.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <io/File.h>
#include <plots/All.h>
#include <settings/All.h>

#include <string>

using namespace ausaxs;

int cli_rigidbody(int argc, char const *argv[]) {
    settings::grid::scaling = 2;
    settings::grid::cubic = true;
    settings::general::verbose = true;
    settings::general::output = "output/rigidbody/";
    bool save_settings = false;

    io::File script, settings;
    CLI::App app{"Perform rigid-body optimization."};
    app.fallthrough();
    auto p_script = app.add_option("input_script", script, "Path to the rigid-body configuration script.")->check(CLI::ExistingFile);
    app.add_flag_callback("--licence",    [] () {console::print_text(constants::licence); exit(0);}, "Print the licence.");
    app.add_flag_callback("-v,--version", [] () {console::print_text(constants::version); exit(0);}, "Print the AUSAXS version.");
    app.add_flag("--allow-unknown-atoms", settings::molecule::allow_unknown_atoms, 
        "Allow processing files with unknown atoms. Use only if you understand the implications.")
        ->default_val(settings::molecule::allow_unknown_atoms);
    app.add_flag("--allow-unknown-residues", settings::molecule::allow_unknown_residues,
        "Allow processing files with unknown residues. Use only if you understand the implications.")
        ->default_val(settings::molecule::allow_unknown_residues);
    app.add_flag("--offline", settings::general::offline, "Run the program in offline mode. This will prevent any network requests.")
        ->default_val(settings::general::offline);
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads);    

    // config subcommands
    auto sub_config = app.add_subcommand("config", "See and set additional options for the configuration.");
    auto p_settings = sub_config->add_option("--file,-f", settings, "The configuration file to use.")->check(CLI::ExistingFile);
    sub_config->add_flag("--save", save_settings, "Save the settings to a file.");
    sub_config->add_flag_callback("--log", [] () {logging::start("saxs_fitter");}, "Enable logging to a file.");

    // data subcommands
    auto sub_data = app.add_subcommand("data", "See and set additional options for the SAXS data.");
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
    sub_data->add_option_function<std::string>("--unit,-u", [] (const std::string& s) {settings::detail::parse_option("unit", {s});}, 
        "The unit of the q values in the measurement file. Options: A, nm.");
    sub_data->add_option("--skip", settings::axes::skip, "Number of points to skip in the measurement file.")->default_val(settings::axes::skip);
    sub_data->add_flag("--rebin", settings::flags::data_rebin, "Rebin the data to increase the information content of each data point.")->default_val(settings::flags::data_rebin);
    sub_data->add_flag("--weighted-bins", settings::hist::weighted_bins, "Decides whether weighted bins are used.")->default_val(settings::hist::weighted_bins);

    // molecule subcommands
    auto sub_mol = app.add_subcommand("molecule", "See and set additional options for the molecular structure file.");
    sub_mol->add_flag("--center,!--no-center", settings::molecule::center, 
        "Decides whether the protein will be centered.")->default_val(settings::molecule::center);
    sub_mol->add_flag("--use-occupancy,!--ignore-occupancy", settings::molecule::use_occupancy, 
        "Decides whether the atomic occupancies from the file will be used.")->default_val(settings::molecule::use_occupancy);

    // hydrogen subcommands
    auto sub_hydrogen = app.add_subcommand("hydrogens", "See and set additional options for the handling of hydration atoms.");
    sub_hydrogen->add_flag("--keep,!--discard", settings::general::keep_hydrogens, "Keep or discard hydrogens from the structure file.")->default_val(settings::general::keep_hydrogens);

    // grid subcommands
    auto sub_grid = app.add_subcommand("grid", "See and set additional options for the grid calculations.");
    sub_grid->add_option("--width,-w", settings::grid::cell_width, 
        "The distance between each grid point in Ångström. Lower widths increase the precision."
    )->default_val(settings::grid::cell_width);

    app.final_callback([&] () {if (save_settings) {
        settings::write("settings.txt");
        console::print_info("Settings saved to settings.txt in current directory.");
    }});

    CLI11_PARSE(app, argc, argv);

    console::print_info("Running AUSAXS " + std::string(constants::version));
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::flags::init_histogram_manager = false;
    logging::start("rigidbody");

    //###################//
    //### PARSE INPUT ###//
    //###################//
    try {
        // rigid-body optimization is only available through the sequencer configuration format
        if (p_script->count() == 0) {
            throw except::missing_option("rigidbody: A configuration script must be supplied.");
        }
        if (!constants::filetypes::config.check(script)) {
            throw except::invalid_argument(
                "rigidbody: The input file \"" + script.str() + "\" is not a configuration script. "
                "Rigid-body optimization must be described in a configuration script."
            );
        }

        // if a settings file was provided, read it and re-parse the command line arguments so they take priority
        if (p_settings->count() != 0) {
            settings::read(settings);
            CLI11_PARSE(app, argc, argv);
        }

        auto res = rigidbody::sequencer::SequenceParser().parse_file(script)->execute();
        fitter::FitReporter::report(res.get());
        fitter::FitReporter::save(res.get(), settings::general::output + "fit.txt");
        res->curves.save(settings::general::output + "ausaxs.fit", "chi2=" + std::to_string(res->fval/res->dof) + " dof=" + std::to_string(res->dof));
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        throw e;
    }
    return 0;
}
