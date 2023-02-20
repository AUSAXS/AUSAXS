#include <CLI/CLI.hpp>

#include <vector>
#include <string>

#include <data/Body.h>
#include <data/Protein.h>
#include <rigidbody/RigidBody.h>
#include <data/BodySplitter.h>
#include <fitter/FitReporter.h>
#include <plots/all.h>

int main(int argc, char const *argv[]) { 
    setting::grid::scaling = 1;
    setting::grid::cubic = true;

    CLI::App app{"Rigid-body optimization."};
    std::string input_structure, input_measurement, output, placement_strategy;
    app.add_option("input_s", input_structure, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", input_measurement, "Path to the measuremed data.")->required()->check(CLI::ExistingFile);
    app.add_option("output", output, "Path to save the hydrated file at.")->required();
    app.add_option("--reduce,-r", setting::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", setting::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width", setting::axes::distance_bin_width, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--placement_strategy", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--qmin", setting::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", setting::axes::qmax, "Upper limit on used q values from measurement file.");
    app.add_flag("--center,!--no-center", setting::protein::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", setting::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
    CLI11_PARSE(app, argc, argv);

    // parse strategy
    if (placement_strategy == "Radial") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::JanStrategy;}

    std::vector<int> splits = {9, 99};
    Protein protein = BodySplitter::split("data/LAR1-2/LAR1-2.pdb", splits);
    RigidBody body(protein);

    protein.save(setting::general::output + "initial.pdb");
    body.optimize(input_measurement);
    protein.save(setting::general::output + "optimized.pdb");
    HydrationFitter fitter(input_measurement, protein.get_histogram());
    auto res = fitter.fit();
    FitReporter::report(res);
    FitReporter::save(res, setting::general::output + "fit.txt");
    plots::PlotIntensityFit::quick_plot(res, setting::general::output + "fit.png");
    return 0;
}