#include <CLI/CLI.hpp>

#include <iostream>
#include <filesystem>

#include <plots/All.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <data/Protein.h>

using std::string;

int main(int argc, char const *argv[]) {
    settings::protein::use_effective_charge = false;
    settings::fit::verbose = true;

    io::ExistingFile mapfile, pdbfile;
    CLI::App app{"Generate a simulated SAXS measurement from an EM map."};
    app.add_option("input-map", mapfile, "Path to the map file. This is only used to generate a suiting output folder.")->required()->check(CLI::ExistingFile);
    app.add_option("input-pdb", pdbfile, "Path to the PDB structure map.")->required()->check(CLI::ExistingFile);
    CLI11_PARSE(app, argc, argv);

    em::ImageStack map(mapfile); 
    Protein protein(pdbfile);
    SimpleDataset msim = protein.simulate_dataset();
    msim.save(io::File(mapfile.directory() + "simulated_SAXS.dat"));
}