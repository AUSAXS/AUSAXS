#include <CLI/CLI.hpp>

#include <grid/Grid.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <settings/All.h>

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);

    io::Folder input;
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("folder", input, "Path to the MD SAXS folder.")->required()->check(CLI::ExistingDirectory);
    CLI11_PARSE(app, argc, argv);

    settings::grid::scaling = 1;
    settings::grid::rvol = 3; // water + atom
    settings::general::output = "output/";

    io::ExistingFile env_unordered(input + "excludedvolume_0.pdb");
    io::ExistingFile env_ordered(input + "prot+solvlayer_0.pdb");

    data::Molecule unordered(env_unordered);
    data::Molecule ordered(env_ordered);
    auto ordered_layer = ordered.get_waters();

    ordered.clear_hydration();
    auto grid = ordered.get_grid();
    std::vector<data::record::Water> random_waters;
    for (auto& water : unordered.get_waters()) {
        auto bin = grid->to_bins(water.get_coordinates());
        if (grid->grid.is_empty(bin.x(), bin.y(), bin.z())) {
            random_waters.emplace_back(water);
        }
    }

    data::Molecule(std::vector<data::record::Atom>{}, ordered_layer).save(settings::general::output + "ordered.pdb");
    data::Molecule(std::vector<data::record::Atom>{}, random_waters).save(settings::general::output + "random.pdb");
    return 0;
}