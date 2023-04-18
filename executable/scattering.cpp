#include <CLI/CLI.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <data/Body.h>
#include <data/Protein.h>
#include <fitter/Fit.h>
#include <plots/all.h>
#include <fitter/FitReporter.h>

int main(int argc, char const *argv[]) {
    CLI::App app{"Calculate the scattering from a pdb structure."};

    std::string pdb, placement_strategy = "Radial";
    bool use_existing_hydration = false, remove_h = false;
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", setting::general::output, "Path to save the generated figures at.")->default_val("output/scattering/");
    app.add_option("--reduce,-r", setting::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,--gw", setting::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width,--bw", setting::axes::distance_bin_width, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--placement_strategy,--ps", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--radius_a,--ra", setting::grid::ra, "Radius of the protein atoms.");
    app.add_option("--radius_h,--rh", setting::grid::rh, "Radius of the hydration atoms.");
    app.add_option("--qlow", setting::axes::qmin, "Lower limit on the q values. Default: 1e-4");
    app.add_option("--qhigh", setting::axes::qmax, "Upper limit on the q values. Default: 0.5");
    app.add_flag("--effective-charge,!--no-effective-charge", setting::protein::use_effective_charge, "Decides whether the effective atomic charge will be used. Default: true.");
    auto opt = app.add_flag("--use-existing-hydration,!--no-use-existing-hydration", use_existing_hydration, "Decides whether the hydration layer will be generated from scratch or if the existing one will be used. Default: false.");
    app.add_flag("--remove-h", remove_h, "Remove all hydrogens from the structure. Default: false.")->excludes(opt);
    CLI11_PARSE(app, argc, argv);

    //####################//
    //### PARSE INPUT ###//
    //####################//
    setting::general::output += utility::stem(pdb) + "/";

    // validate input
    if (!constants::filetypes::structure.validate(pdb)) {
        throw except::invalid_argument("Unknown PDB extensions: " + pdb);
    }

    // parse strategy
    if (placement_strategy == "Radial") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::JanStrategy;}

    //######################//
    //### ACTUAL PROGRAM ###//
    //######################//
    Protein protein(pdb);
    if (!use_existing_hydration || protein.hydration_atoms.empty()) {
        protein.generate_new_hydration();
    }
    hist::ScatteringHistogram h = protein.get_histogram();
    plots::PlotDistance::quick_plot(h, setting::general::output + "p(r)." + setting::plot::format);
    plots::PlotIntensity::quick_plot(h, setting::general::output + "scattering." + setting::plot::format);

    std::cout << "I(" << h.q[0] << ") = " << h.p[0] << std::endl;

    return 0;
}