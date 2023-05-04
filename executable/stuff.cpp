#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/Protein.h>
#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <fitter/FitReporter.h>
#include <plots/all.h>
#include <settings/All.h>
#include <io/File.h>
#include <rigidbody/sequencer/Sequencer.h>

int main(int argc, char const *argv[]) {
    settings::grid::scaling = 2;
    settings::grid::cubic = true;
    settings::general::verbose = true;
    settings::axes::distance_bin_width = 0.1;

    CLI::App app{"Rigid-body optimization."};
    io::File pdb, mfile, settings;
    std::vector<unsigned int> constraints;
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", mfile, "Path to the measuremed data.")->required()->check(CLI::ExistingFile);
    app.add_option("output", settings::general::output, "Path to save the hydrated file at.")->default_val("output/rigidbody/");
    auto p_cal = app.add_option("--calibrate", settings::rigidbody::detail::calibration_file, "Path to the calibration data.")->expected(0, 1)->check(CLI::ExistingFile);
    app.add_option("--reduce,-r", settings::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", settings::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width", settings::axes::distance_bin_width, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--iterations", settings::rigidbody::iterations, "Maximum number of iterations. Default: 1000.");
    app.add_option("--constraints", settings::rigidbody::detail::constraints, "Constraints to apply to the rigid body.");
    app.add_flag("--center,!--no-center", settings::protein::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", settings::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
    CLI11_PARSE(app, argc, argv);
    
    //###################//
    //### PARSE INPUT ###//
    //###################//
    settings::general::output += mfile.stem() + "/";

    // if a settings file was provided
    if (p_settings->count() != 0) {
        settings::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (settings::discover(std::filesystem::path(mfile).parent_path().string())) {
            CLI11_PARSE(app, argc, argv);
        }
    }
    if (settings::rigidbody::detail::constraints.empty()) {
        throw except::missing_option("rigidbody: Constraints must be supplied. Use --constraints to specify them.");
    }

    rigidbody::sequencer::Sequencer(mfile, rigidbody::BodySplitter::split(pdb, settings::rigidbody::detail::constraints))
        .body_select_strategy(settings::rigidbody::BodySelectStrategyChoice::RandomConstraintSelect)
        .parameter_strategy(settings::rigidbody::ParameterGenerationStrategyChoice::RotationsOnly)
            .decay_strategy(settings::rigidbody::DecayStrategyChoice::Exponential)
            .amplitude(0.5)
        .transform_strategy(settings::rigidbody::TransformationStrategyChoice::RigidTransform)
        .loop(5)
    .execute();

    return 0;
}