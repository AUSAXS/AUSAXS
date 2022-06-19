#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::fit::q_high = 2;

    // string mapfile = "data/lysozyme/emd_23957.map";
    // string mapfile = "data/A2M/emd_12747.map";
    // string mapfile = "data/A2M/emd_12748.map";
    // string mapfile = "data/A2M/emd_12752.map"; // tryp
    // string mapfile = "data/A2M/emd_12753.map"; // tryp
    // string mapfile = "data/SHOC2/emd_25044.map";
    // string mapfile = "data/SHOC2/emd_26667.map";
    string mapfile = "data/flipped_ns_igefceria.mrc";

    // string pdbfile = "data/2epe.pdb";
    string pdbfile = "data/native.pdb";

    // string mfile = "data/2epe.RSR";
    // string mfile = "data/A2M_native.RSR";
    // string mfile = "data/A2M_tryp.RSR";
    string mfile = "data/shoc2.dat";

    em::ImageStack map(mapfile); 

    //* STRUCTURE FIT
    if (false) {
        string path = "figures/em/structure_fit/" + utility::stem(mapfile) + "/" + utility::stem(pdbfile) + "/";

        Protein pdb(pdbfile);
        auto pdb_h = pdb.get_histogram();
        auto res = map.fit(pdb_h);
        FitReporter::report(res);
        FitReporter::save(path + "report.txt", res);

        plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
        plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

        auto scan = map.cutoff_scan(100, pdb_h);
        plots::PlotDataset::quick_plot(scan, path + "scan.pdf");
    }

    //* MEASUREMENT FIT
    if (false) {
        string path = "figures/em/measurement_fit/" + utility::stem(mfile) + "/" + utility::stem(mapfile) + "/";

        auto res = map.fit(mfile);
        FitReporter::report(res);
        FitReporter::save(path + "report.txt", res);

        plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
        plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

        auto scan = map.cutoff_scan(100, mfile);
        plots::PlotDataset::quick_plot(scan, path + "scan.pdf");
    }

    //* GENERATE INTENSITY & PDB
    if (true) {
        // unsigned int c;
        // for(const em::Image& image : map.images()) {
        //     plots::PlotImage::quick_plot(image, "figures/em/images/" + utility::stem(mapfile) + "/" + std::to_string(c++) + ".png");
        // }
        plots::PlotIntensity::quick_plot(map.get_histogram(map.level(3)), "figures/em/cut/intensity.pdf");
        map.save("data/output/" + utility::stem(mapfile) + ".pdb", map.level(1));
    }

    return 0;
}