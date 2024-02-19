#include <CLI/CLI.hpp>

#include <string>

#include <data/Body.h>
#include <data/Molecule.h>
#include <settings/All.h>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::string s_pdb, s_mfile, s_settings, placement_strategy = "radial", histogram_manager = "hmmt"; // not using partial histograms has a slightly smaller overhead

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input", s_pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/intensity_fitter/")->group("General options");
    app.add_flag("--center,!--no-center", settings::molecule::center, "Decides whether the protein will be centered.")->default_val(settings::molecule::center)->group("Protein options");

    app.add_option("--reduce,-r", settings::grid::water_scaling, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.")->default_val(settings::grid::water_scaling)->group("Advanced options");
    app.add_option("--grid_width,--gw", settings::grid::width, "The distance between each grid point in Ångström. Lower widths increase the precision.")->default_val(settings::grid::width)->group("Advanced options");
    app.add_option("--placement_strategy,--ps", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.")->default_val(placement_strategy)->group("Advanced options");
    app.add_option("--exv_radius,--er", settings::grid::exv_radius, "The radius of the excluded volume sphere used for the grid-based excluded volume calculations in Ångström.")->default_val(settings::grid::exv_radius)->group("Advanced options");
    CLI11_PARSE(app, argc, argv);

    // parse strategy
    if (placement_strategy == "Radial") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::JanStrategy;}

    data::Molecule protein(s_pdb);
    protein.generate_new_hydration();
    protein.save(settings::general::output + "hydration.pdb");
    return 0;
}