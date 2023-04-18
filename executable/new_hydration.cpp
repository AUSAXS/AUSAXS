#include <CLI/CLI.hpp>

#include <vector>
#include <string>

#include <data/Body.h>
#include <data/Protein.h>
#include <utility/Settings.h>

int main(int argc, char const *argv[]) {
    CLI::App app{"Generate a new hydration layer for a given input data file."};
    setting::grid::scaling = 2;

    std::string input, output, placement_strategy;
    app.add_option("input", input, "Path to the data file.")->required()->check(CLI::ExistingFile);
    app.add_option("output", output, "Path to save the hydrated file at.")->required();
    app.add_option("--reduce,-r", setting::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", setting::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--placement_strategy", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--radius_a", setting::grid::ra, "Radius of the protein atoms.");
    app.add_option("--radius_h", setting::grid::rh, "Radius of the hydration atoms.");
    app.add_flag("--center,!--no-center", setting::protein::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", setting::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
    CLI11_PARSE(app, argc, argv);

    // parse strategy
    if (placement_strategy == "Radial") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::JanStrategy;}

    Protein protein(input);
    protein.generate_new_hydration();
    protein.save(output);
    return 0;
}