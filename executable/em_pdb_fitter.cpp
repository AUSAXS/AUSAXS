#include <iostream>
#include <filesystem>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::fit::verbose = true;
    setting::em::sample_frequency = 1;
    setting::em::alpha_levels = {1, 10};

    // check that we have at least one argument
    if (argc == 5) {
        setting::em::alpha_levels = {std::stod(argv[3]), std::stod(argv[4])};
    } else if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "<mapfile> <pdbfile>" << std::endl;
        return 1;
    }

    // load the input files
    string mapfile = argv[1];
    string pdbfile = argv[2];
    setting::plot::path = "figures/em_pdb_fitter/" + utility::stem(mapfile) + "/";
    string mfile = setting::plot::path + "temp.dat";

    // load the map and protein
    em::ImageStack map(mapfile); 
    Protein protein(pdbfile);

    // simulate a SAXS measurement
    SimpleDataset msim = protein.simulate_dataset(false);
    msim.save(mfile);
    auto res = map.fit(mfile);
    std::filesystem::remove(mfile);

    plots::PlotIntensityFit::quick_plot(res, setting::plot::path + "fit." + setting::plot::format);
    plots::PlotIntensityFitResiduals::quick_plot(res, setting::plot::path + "residuals." + setting::plot::format);

    std::cout << "DOF: " << res->dof << std::endl;
    FitReporter::report(res);
    FitReporter::save(res, setting::plot::path + "report.txt");
}