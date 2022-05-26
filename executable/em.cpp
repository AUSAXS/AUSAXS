#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    string mapfile = "data/A2M_ma.map";
    string pdbfile = "data/A2M_ma.pdb";
    em::ImageStack map(mapfile); 
    std::cout << map.get_header()->to_string() << std::endl;

    plots::PlotIntensity::quick_plot(map.get_histogram(2), "figures/em/intensity.pdf");

    Protein pdb(pdbfile);
    auto res = map.fit(pdb.get_histogram());
    FitReporter::report(res);
    return 0;
}