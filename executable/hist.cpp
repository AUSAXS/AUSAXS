#include <CLI/CLI.hpp>

#include <vector>
#include <string>

#include <data/Body.h>
#include <data/Protein.h>
#include <plots/PlotDistance.h>
#include <plots/PlotIntensity.h>
#include <settings/All.h>

int main(int argc, char const *argv[]) { 
    CLI::App app{"Generate a distance histogram and a scattering intensity plot for a given input data file."};
    app.prefix_command(false);

    settings::grid::scaling = 2;
    settings::axes::distance_bin_width = 0.1;

    std::string input, output, placement_strategy;
    app.add_option("input", input, "Path to the data file.")->required()->check(CLI::ExistingFile);
    app.add_option("output", output, "Path to save the hydrated file at.")->required();
    app.add_option("--reduce,-r", settings::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", settings::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width", settings::axes::distance_bin_width, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--placement_strategy", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--radius_a", settings::grid::ra, "Radius of the protein atoms.");
    app.add_option("--radius_h", settings::grid::rh, "Radius of the hydration atoms.");
    app.add_flag("--center,!--no-center", settings::protein::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", settings::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
    CLI11_PARSE(app, argc, argv);

    // parse strategy
    if (placement_strategy == "Radial") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::JanStrategy;}

    settings::protein::use_effective_charge = false;

    settings::general::keep_hydrogens = true;
    Protein protein(input);
    protein.clear_hydration();
    protein.save("tmp/TEST.PDB");
    protein.generate_new_hydration();
    hist::ScatteringHistogram d = protein.get_histogram();

    // Distance plot
    plots::PlotDistance d_plot(d, output + "distance." + settings::plots::format);

    // Debye scattering intensity plot
    // plots::PlotIntensity i_plot(d);
    // i_plot.plot_guinier_approx(d);
    // i_plot.save(output + "intensity." + settings::plot::format);

    std::cout << "Protein size: " << protein.atom_size() << std::endl;

    double dmax = 0;
    auto atoms = protein.atoms();
    for (auto &atom1 : atoms) {
        for (auto &atom2 : atoms) {
            if (atom1 != atom2) {
                double d = atom1.distance(atom2);
                if (d > dmax) {
                    dmax = d;
                }
            }
        }
    }
    std::cout << "Maximum distance: " << dmax << std::endl;
    return 0;
}