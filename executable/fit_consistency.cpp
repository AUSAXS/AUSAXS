#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::grid::psc = setting::grid::PlacementStrategyChoice::RadialStrategy;
    setting::em::sample_frequency = 2;

    // check that we have at least one argument
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << "<mapfile> <pdbfile> <mfile>" << std::endl;
        return 1;
    }

    // load the input files
    string mapfile = argv[1];
    string pdbfile = argv[2];
    string mfile = argv[3];
    string path = "figures/consistency/" + utility::stem(mapfile) + "/" + utility::stem(pdbfile) + "/";
    em::ImageStack map(mapfile); 
    Protein protein(pdbfile);

    // fit the measurement to the protein
    protein.generate_new_hydration();
    auto res1 = protein.fit(mfile);
    plots::PlotIntensityFit::quick_plot(res1, path + "intensity_fit." + setting::figures::format);
    plots::PlotIntensityFitResiduals::quick_plot(res1, path + "residuals." + setting::figures::format);
    FitReporter::report(res1);

    // fit the measurement to the map
    auto res2 = map.fit(mfile);
    plots::PlotIntensityFit::quick_plot(res2, path + "intensity_fit_map." + setting::figures::format);
    plots::PlotIntensityFitResiduals::quick_plot(res2, path + "residuals_map." + setting::figures::format);
    FitReporter::report(res2);
}