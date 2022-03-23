// includes
#include <vector>
#include <string>

#include <data/Body.h>
#include <data/Protein.h>
#include <plots/PlotDistance.h>
#include <plots/PlotIntensity.h>
#include <CLI11.hpp>

int main(int argc, char const *argv[]) { 
    CLI::App app{"Generate a distance histogram and a scattering intensity plot for a given input data file."};

    std::string input, output, placement_strategy;
    app.add_option("input", input, "Path to the data file.")->required()->check(CLI::ExistingFile);
    app.add_option("output", output, "Path to save the hydrated file at.")->required();
    app.add_option("--reduce,-r", setting::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,-w", setting::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width", setting::axes::scattering_intensity_plot_binned_width, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--placement_strategy", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--radius_a", setting::grid::ra, "Radius of the protein atoms.");
    app.add_option("--radius_h", setting::grid::rh, "Radius of the hydration atoms.");
    app.add_flag("--center,!--no-center", setting::protein::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", setting::protein::use_effective_charge, "Decides whether the protein will be centered. Default: true.");
    CLI11_PARSE(app, argc, argv);

    // parse strategy
    if (placement_strategy == "Radial") {setting::grid::psc = setting::grid::PlacementStrategyChoice::RadialStrategy;}
    else if (placement_strategy == "Axes") {setting::grid::psc = setting::grid::PlacementStrategyChoice::AxesStrategy;}
    else if (placement_strategy == "Jan") {setting::grid::psc = setting::grid::PlacementStrategyChoice::JanStrategy;}

    setting::protein::use_effective_charge = false;
    setting::axes::scattering_intensity_plot_axis = {1000, 0.001, 100.001};
    setting::axes::scattering_intensity_plot_binned_width = 0.01;

    Protein protein(input);
    protein.generate_new_hydration();
    ScatteringHistogram d = protein.get_histogram();

    // Distance plot
    plots::PlotDistance d_plot(d);
    d_plot.save(output + "distances." + setting::figures::format); 

    // Debye scattering intensity plot
    plots::PlotIntensity i_plot(d);
    i_plot.plot_guinier_approx();
    i_plot.save(output + "intensity." + setting::figures::format);

    std::cout << "Protein size: " << protein.atom_size() << std::endl;
    return 0;
}