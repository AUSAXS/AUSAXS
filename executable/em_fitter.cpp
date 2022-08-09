#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "<mapfile> <mfile>" << std::endl;
        return 1;
    }

    string mapfile = argv[1];
    string mfile = argv[2];

    std::cout << "Performing EM fit with map " << mapfile << " and measurement " << mfile << std::endl;

    em::ImageStack map(mapfile); 

    // Fit the measurements to the EM density map.

    string path = "figures/em_fitter/" + utility::stem(mfile) + "/" + utility::stem(mapfile) + "/";

    auto res = map.fit(mfile);
    FitReporter::report(res);
    FitReporter::save(path + "report.txt", res);

    plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
    plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

    // auto scan = map.cutoff_scan(100, mfile);
    // plots::PlotDataset::quick_plot(scan, path + "scan.pdf");
    return 0;
}