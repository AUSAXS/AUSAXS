// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_rigidbody.h>
#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/DefaultOptimizer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <fitter/FitReporter.h>
#include <fitter/SmartFitter.h>
#include <constants/Constants.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <io/File.h>
#include <plots/All.h>
#include <settings/All.h>

#include <vector>
#include <string>

using namespace ausaxs;

void add_calibration(rigidbody::Rigidbody& rigidbody, const io::ExistingFile& mfile);
int cli_rigidbody(int argc, char const *argv[]) {
    settings::grid::scaling = 2;
    settings::grid::cubic = true;
    settings::general::verbose = true;
    bool save_settings = false;

    io::File pdb, mfile, settings, config;
    std::vector<unsigned int> constraints;
    CLI::App app{"[EXPERIMENTAL] Perform rigid-body optimization."};
    app.fallthrough();
    auto input_s = app.add_option("input_structure", pdb, "Path to the structure file.")->check(CLI::ExistingFile);
    auto input_m = app.add_option("input_measurement", mfile, "Path to the measured SAXS data.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Output folder to write the results to.")->default_val("output/rigidbody/")->group("General options");
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

    // hidden options group
    app.add_flag("--weighted-bins", settings::hist::weighted_bins, 
        "Decides whether the weighted bins will be used."
    )->default_val(settings::hist::weighted_bins)->group("");

    auto p_cal = app.add_option("--calibrate", settings::rigidbody::detail::calibration_file, "Path to the calibration data.")->expected(0, 1)->check(CLI::ExistingFile)->group("");
    app.add_option("--constraints", settings::rigidbody::detail::constraints, "Constraints to apply to the rigid body.")->group("");
    CLI11_PARSE(app, argc, argv);

    console::print_info("Running AUSAXS " + std::string(constants::version));
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::flags::init_histogram_manager = false;
    logging::start("rigidbody");

    //###################//
    //### PARSE INPUT ###//
    //###################//
    try {
        settings::general::output += mfile.stem() + "/";

        // check if pdb is a config script
        if (constants::filetypes::config.check(pdb)) {
            auto res = rigidbody::sequencer::SequenceParser().parse(pdb)->execute();
            fitter::FitReporter::report(res.get());
            fitter::FitReporter::save(res.get(), settings::general::output + "fit.txt");
            res->curves.save(settings::general::output + "ausaxs.fit", "chi2=" + std::to_string(res->fval/res->dof) + " dof=" + std::to_string(res->dof));
            return 0;
        }

        // if a settings file was provided
        if (p_settings->count() != 0) {
            settings::read(settings);        // read it
            CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
        } else {                            // otherwise check if there is a settings file in the same directory
            if (settings::discover(std::filesystem::path(mfile.path()).parent_path().string())) {
                CLI11_PARSE(app, argc, argv);
            }
        }
        if (settings::rigidbody::detail::constraints.empty()) {
            throw except::missing_option("Constraints must be supplied. Use --constraints to specify them.");
        }
        settings::validate_settings();

        rigidbody::Rigidbody rigidbody = rigidbody::BodySplitter::split(pdb, settings::rigidbody::detail::constraints);
        if (p_cal->count() != 0) { // calibration file was provided
            add_calibration(rigidbody, mfile);
        }
        auto res = rigidbody::default_optimize(&rigidbody, mfile);
        fitter::FitReporter::report(res.get());
        fitter::FitReporter::save(res.get(), settings::general::output + "fit.txt", argc, argv);
        res->curves.save(settings::general::output + "ausaxs.fit", "chi2=" + std::to_string(res->fval/res->dof) + " dof=" + std::to_string(res->dof));
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        throw e;
    }
    return 0;
}

void add_calibration(rigidbody::Rigidbody& rigidbody, const io::ExistingFile& mfile) {
    if (settings::rigidbody::detail::calibration_file.empty()) {
        // look for calibration file in the same directory as the measurement file
        for (const auto& f : mfile.directory().files()) {
            if (constants::filetypes::saxs_data.check(f)) {
                settings::rigidbody::detail::calibration_file = f.path();
                std::cout << "\tUsing calibration file: \"" << f.path() << "\"" << std::endl;
                break;
            }
        }
    }
    if (settings::rigidbody::detail::calibration_file.empty()) {
        throw except::missing_option("rigidbody: Default calibration file not found. Use --calibrate <path> to specify it.");
    }

    settings::general::output += "calibrated/";
    rigidbody.molecule.generate_new_hydration();
    fitter::SmartFitter fitter({settings::rigidbody::detail::calibration_file}, rigidbody.molecule.get_histogram());
    auto res = fitter.fit();
    if (settings::general::verbose) {
        std::cout << "Calibration results:" << std::endl;
        fitter::FitReporter::report(res.get());
    }
    throw std::runtime_error("rigidbody: Calibration is currently disabled.");
    // rigidbody.apply_calibration(std::move(res));
}