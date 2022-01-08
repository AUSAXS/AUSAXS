// includes
#include <vector>
#include <string>

#include "data/Body.h"
#include "data/Protein.h"
#include "plot_style.h"
#include "plots/PlotDistance.h"
#include "plots/PlotIntensity.h"
#include "CLI11.hpp"

int main(int argc, char const *argv[]) { 
    CLI::App app{"Generate a distance histogram and a scattering intensity plot for a given input data file."};

    std::string input, output, placement_strategy;
    app.add_option("input", input, "Path to the data file.")->required()->check(CLI::ExistingFile);
    app.add_option("output", output, "Path to save the hydrated file at.")->required();
    app.add_option("--reduce,-r", setting::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--width,-w", setting::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--placement_strategy", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--radius_a", setting::grid::ra, "Radius of the protein atoms.");
    app.add_option("--radius_h", setting::grid::rh, "Radius of the hydration atoms.");
    CLI11_PARSE(app, argc, argv);

    // setting::axes::scattering_intensity_plot_binned_width = 0.5;
    // setting::figures::format = "png";

    Protein protein(input);
    protein.generate_new_hydration();
    shared_ptr<ScatteringHistogram> d = protein.get_histogram();

    // Distance plot
    PlotDistance d_plot(d);
    d_plot.save(output + "distances." + setting::figures::format); 

    // Debye scattering intensity plot
    PlotIntensity i_plot(d);
    i_plot.save(output + "intensity." + setting::figures::format);
    return 0;
}