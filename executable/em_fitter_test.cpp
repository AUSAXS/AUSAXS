#include <iostream>

#include <plots/All.h>
#include <data/Molecule.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>

using std::string;

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    settings::em::mass_axis = true;
    settings::em::hydrate = true;
    settings::fit::verbose = true;
    settings::em::alpha_levels = {0.5, 5};
    settings::hist::weighted_bins = true;
    settings::molecule::use_effective_charge = true;
    settings::molecule::implicit_hydrogens = true;
    settings::em::fixed_weights = true;
    settings::general::threads /= 2;

    // check that we have all three arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "<mapfile> <pdbfile>" << std::endl;
        return 1;
    }

    // load the input files
    io::ExistingFile mapfile = argv[1];
    io::ExistingFile pdbfile = argv[2];
    settings::general::output += "em_fitter_test/" + mapfile.stem() + "/" + pdbfile.stem() + "/";

    // load the map and protein
    data::Molecule protein(pdbfile);
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    em::ImageStack map(mapfile); 

    // simulate saxs data from protein
    auto saxs_path = settings::general::output + "simulated_saxs.dat";
    {
        auto saxs = protein.simulate_dataset();
        saxs.save(saxs_path);
    }

    // fit the map to the data
    auto fit = map.fit(saxs_path);
    fitter::FitReporter::report(fit.get());
    fitter::FitReporter::save(fit.get(), settings::general::output + "report.txt");
}