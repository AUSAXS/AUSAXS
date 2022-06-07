#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;

    string mapfile = "data/A2M_ma.map";
    string pdbfile = "data/A2M_ma.pdb";
    em::ImageStack map(mapfile); 
    std::cout << map.get_header()->to_string() << std::endl;

    plots::PlotIntensity::quick_plot(map.get_histogram(0.05), "figures/em/intensity.pdf");

    Protein pdb(pdbfile);
    auto pdb_h = pdb.get_histogram();
    auto res = map.fit(pdb_h, mini::Parameter("cutoff", 0.05, {0.01, 0.1}));
    FitReporter::report(res);

    auto scan = map.cutoff_scan(Axis(10, 0, 0.1), pdb_h);
    plots::PlotDataset::quick_plot(scan, "figures/em/scan.pdf");
    return 0;
}