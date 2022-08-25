#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = true;
    setting::fit::verbose = false;
    setting::grid::placement_strategy = setting::grid::PlacementStrategy::RadialStrategy;
    setting::em::sample_frequency = 2;

    // check that we have at least one argument
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "<mapfile> <pdbfile>" << std::endl;
        return 1;
    }

    // load the input files
    string mapfile = argv[1];
    string pdbfile = argv[2];
    string path = "figures/fit_consistency2/" + utility::stem(mapfile) + "/" + utility::stem(pdbfile) + "/";
    string mfile = path + "temp.dat";

    // load the map and protein
    em::ImageStack map(mapfile); 
    Protein protein(pdbfile);

    unsigned int evals = 100;
    Dataset data({"dof", "chi2_struct", "chi2_map", "cutoff"});
    for (unsigned int i = 0; i < evals; i++) {
        std::cout << "Starting iteration " << i+1 << " of " << evals << std::endl;

        // simulate a SAXS measurement
        SimpleDataset msim = protein.simulate_dataset();
        msim.save(mfile);

        // fit the measurement to the protein
        auto res1 = protein.fit(mfile);

        // fit the measurement to the map
        auto res2 = map.fit(mfile);

        // add to output file
        data.push_back({double(res1->dof), res1->fval, res2->fval, res2->get_parameter("cutoff")});
    }
    data.save(path + "out.txt");
}