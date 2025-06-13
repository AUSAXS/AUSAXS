// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/BodySplitter.h>
#include <fitter/FitReporter.h>
#include <fitter/SmartFitter.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <constants/Constants.h>
#include <utility/Console.h>
#include <io/File.h>
#include <plots/All.h>
#include <settings/All.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/DefaultOptimizer.h>

#include <vector>
#include <string>

using namespace ausaxs;

void add_calibration(rigidbody::RigidBody& rigidbody, const io::ExistingFile& mfile);
int main(int argc, char const *argv[]) { 
    settings::grid::scaling = 2;
    settings::grid::cubic = true;
    settings::grid::min_bins = 1000;
    settings::general::verbose = true;

    CLI::App app{"Rigid-body optimization."};
    io::File pdb, mfile, settings, config;
    std::vector<unsigned int> constraints;
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", mfile, "Path to the measuremed data.")->check(CLI::ExistingFile);
    app.add_option("output", settings::general::output, "Path to save the hydrated file at.")->default_val("output/rigidbody/");
    auto p_cal = app.add_option("--calibrate", settings::rigidbody::detail::calibration_file, "Path to the calibration data.")->expected(0, 1)->check(CLI::ExistingFile);
    app.add_option("--reduce,-r", settings::grid::water_scaling, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", settings::grid::cell_width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option_function<std::string>("--placement-strategy,--ps", [] (const std::string& s) {settings::detail::parse_option("placement_strategy", {s});}, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--iterations", settings::rigidbody::iterations, "Maximum number of iterations. Default: 1000.");
    app.add_option("--constraints", settings::rigidbody::detail::constraints, "Constraints to apply to the rigid body.");
    app.add_flag("--center,!--no-center", settings::molecule::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--quit-on-unknown-atom,!--no-quit-on-unknown-atom", settings::molecule::throw_on_unknown_atom, "Decides whether the program will quit if an unknown atom is found. Default: true.");
    CLI11_PARSE(app, argc, argv);

    console::print_info("Running AUSAXS " + std::string(constants::version));

    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::molecule::implicit_hydrogens = false;

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
            throw except::missing_option("rigidbody: Constraints must be supplied. Use --constraints to specify them.");
        }
        settings::validate_settings();

        rigidbody::RigidBody rigidbody = rigidbody::BodySplitter::split(pdb, settings::rigidbody::detail::constraints);
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

void add_calibration(rigidbody::RigidBody& rigidbody, const io::ExistingFile& mfile) {
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
    rigidbody.generate_new_hydration();
    fitter::SmartFitter fitter({settings::rigidbody::detail::calibration_file}, rigidbody.get_histogram());
    auto res = fitter.fit();
    if (settings::general::verbose) {
        std::cout << "Calibration results:" << std::endl;
        fitter::FitReporter::report(res.get());
    }
    rigidbody.apply_calibration(std::move(res));
    // plots::PlotIntensityFit::quick_plot(res.get(), settings::general::output + "calibration.png");
}