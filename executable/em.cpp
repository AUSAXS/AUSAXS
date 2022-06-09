#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;

    // string mapfile = "data/lysozyme/emd_23957.map";
    // mini::Parameter param("cutoff", 4, {1, 10}); // 23957
    string mapfile = "data/A2M/emd_12747.map";
    mini::Parameter param("cutoff", 0.05, {0.005, 0.1}); // 12747
    // string mapfile = "data/A2M/emd_12748.map";
    // mini::Parameter param("cutoff", 0.05, {0.01, 0.06}); // 12748
    em::ImageStack map(mapfile); 

    //* STRUCTURE FIT
    if (false) {
        // string pdbfile = "data/2epe.pdb";
        string pdbfile = "data/native.pdb";
        string path = "figures/em/structure_fit/" + utility::stem(mapfile) + "/" + utility::stem(pdbfile) + "/";

        Protein pdb(pdbfile);
        auto pdb_h = pdb.get_histogram();
        auto res = map.fit(pdb_h, param);
        FitReporter::report(res);
        FitReporter::save(path + "report.txt", res);

        plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
        plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

        auto scan = map.cutoff_scan(Axis(100, param.bounds->min, param.bounds->max), pdb_h);
        plots::PlotDataset::quick_plot(scan, path + "scan.pdf");
    }

    //* MEASUREMENT FIT
    if (true) {
        // string mfile = "data/2epe.RSR";
        string mfile = "data/A2M_native.RSR";
        string path = "figures/em/measurement_fit/" + utility::stem(mfile) + "/" + utility::stem(mapfile) + "/";

        auto res = map.fit(mfile, param);
        FitReporter::report(res);
        FitReporter::save(path + "report.txt", res);

        plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
        plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

        auto scan = map.cutoff_scan(Axis(100, param.bounds->min, param.bounds->max), mfile);
        plots::PlotDataset::quick_plot(scan, path + "scan.pdf");
    }

    return 0;
}