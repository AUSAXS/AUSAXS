#include <CLI/CLI.hpp>

#include <vector>
#include <string>

#include <data/Body.h>
#include <data/Protein.h>
#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <fitter/FitReporter.h>
#include <plots/all.h>
#include <utility/AllSettings.h>

int main(int argc, char const *argv[]) { 
    settings::grid::scaling = 2;
    settings::grid::cubic = true;
    settings::general::verbose = true;

    // settings::rigidbody::pgsc = settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly;
    settings::axes::distance_bin_width = 0.1;

    CLI::App app{"Rigid-body optimization."};
    std::string pdb, mfile, placement_strategy, settings;
    std::vector<unsigned int> constraints;
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", mfile, "Path to the measuremed data.")->required()->check(CLI::ExistingFile);
    app.add_option("output", settings::general::output.value, "Path to save the hydrated file at.")->default_val("output/rigidbody/");
    auto p_cal = app.add_option("--calibrate", settings::rigidbody::detail::calibration_file.value, "Path to the calibration data.")->expected(0, 1)->check(CLI::ExistingFile);
    app.add_option("--reduce,-r", settings::grid::percent_water.value, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", settings::grid::width.value, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width", settings::axes::distance_bin_width.value, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--placement_strategy", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--qmin", settings::axes::qmin.value, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", settings::axes::qmax.value, "Upper limit on used q values from measurement file.");
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--iterations", settings::rigidbody::iterations.value, "Maximum number of iterations. Default: 1000.");
    app.add_option("--constraints", settings::rigidbody::detail::constraints.value, "Constraints to apply to the rigid body.");
    app.add_flag("--center,!--no-center", settings::protein::center.value, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", settings::protein::use_effective_charge.value, "Decides whether the protein will be centered. Default: true.");
    CLI11_PARSE(app, argc, argv);
    
    //###################//
    //### PARSE INPUT ###//
    //###################//
    settings::general::output.value += utility::stem(mfile) + "/";

    // if a settings file was provided
    if (p_settings->count() != 0) {
        settings::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (settings::discover(std::filesystem::path(mfile).parent_path().string())) {
            CLI11_PARSE(app, argc, argv);
        }
    }
    if (settings::rigidbody::detail::constraints.value.empty()) {
        throw except::missing_option("rigidbody: Constraints must be supplied. Use --constraints to specify them.");
    }

    // parse strategy
    if (placement_strategy == "Radial") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::JanStrategy;}
    rigidbody::RigidBody rigidbody = rigidbody::BodySplitter::split(pdb, settings::rigidbody::detail::constraints);

    if (p_cal->count() != 0) {
        if (settings::rigidbody::detail::calibration_file.value.empty()) {
            // look for calibration file in the same directory as the measurement file
            std::vector<std::string> valid_names = {"calibration", "gromacs", "waxs_final"};
            for (const auto& name : valid_names) {
                for (const auto& ext : constants::filetypes::saxs_data.extensions) {
                    std::string path = std::filesystem::path(mfile).parent_path().string() + "/" + name + "." + ext;
                    if (std::filesystem::exists(path)) {
                        settings::rigidbody::detail::calibration_file = path;
                        std::cout << "\tUsing calibration file: \"" << path << "\"" << std::endl;
                        break;
                    }
                }
            }
        }
        if (settings::rigidbody::detail::calibration_file.value.empty()) {
            throw except::missing_option("rigidbody: Default calibration file not found. Use --calibrate <path> to specify it.");
        }

        settings::general::output.value += "calibrated/";
        rigidbody.generate_new_hydration();
        fitter::HydrationFitter fitter(settings::rigidbody::detail::calibration_file, rigidbody.get_histogram());
        auto res = fitter.fit();
        if (settings::general::verbose) {
            std::cout << "Calibration results:" << std::endl;
            fitter::FitReporter::report(res);
        }
        rigidbody.apply_calibration(res);
        plots::PlotIntensityFit::quick_plot(res, settings::general::output + "calibration.png");
    } else {
        settings::general::output.value += "uncalibrated/";
    }

    rigidbody.save(settings::general::output + "initial.pdb");
    rigidbody.optimize(mfile);    
    rigidbody.save(settings::general::output + "optimized.pdb");

    std::shared_ptr<fitter::Fit> res;
    if (!settings::rigidbody::detail::calibration_file.value.empty()) {
        std::shared_ptr<fitter::LinearFitter> fitter = std::make_shared<fitter::LinearFitter>(mfile);
        rigidbody.update_fitter(fitter);
        res = fitter->fit();
    } else {
        std::shared_ptr<fitter::HydrationFitter> fitter = std::make_shared<fitter::HydrationFitter>(mfile, rigidbody.get_histogram());
        res = fitter->fit();
    }
    fitter::FitReporter::report(res);
    fitter::FitReporter::save(res, settings::general::output + "fit.txt");
    plots::PlotIntensityFit::quick_plot(res, settings::general::output + "fit.png");
    return 0;
}