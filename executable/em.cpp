#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    // setting::fit::q_high = 2;
    // setting::axes::scattering_intensity_plot_axis = {1000, 1e-3, 2};
    // setting::axes::scattering_intensity_plot_binned_width = 0.01;

    // string mapfile = "data/lysozyme/emd_23957.map";
    // string mapfile = "data/A2M/emd_12747.map";
    // string mapfile = "data/A2M/emd_12748.map";
    // string mapfile = "data/A2M/emd_12752.map"; // tryp
    // string mapfile = "data/A2M/emd_12753.map"; // tryp
    std::string mapfile = "data/SHOC2/emd_25044.map";
    // string mapfile = "data/SHOC2/emd_26667.map";
    // string mapfile = "data/flipped_ns_igefceria.mrc";
    // string mapfile = "sim/2epe_10.ccp4";

    // string pdbfile = "data/lysozyme/2epe.pdb";
    std::string pdbfile = "data/SHOC2/7sd0.pdb";
    // string pdbfile = "data/native.pdb";

    // string mfile = "data/lysozyme/2epe.RSR";
    // string mfile = "data/A2M_native.RSR";
    // string mfile = "data/A2M_tryp.RSR";
    std::string mfile = "data/SHOC2/7sd0.dat";

    // override the default settings if we have arguments
    if (argc == 4) {
        mapfile = argv[1];
        pdbfile = argv[2];
        mfile = argv[3];
    }
    std::cout << "Performing EM fit with map " << mapfile << " and protein " << pdbfile << " and measurement " << mfile << std::endl;

    em::ImageStack map(mapfile); 

    //* STRUCTURE FIT
    // Fit the entire structure file to the EM density map.
    if (false) {
        std::string path = "figures/em/structure_fit/" + utility::stem(mapfile) + "/" + utility::stem(pdbfile) + "/";

        Protein pdb(pdbfile);
        auto pdb_h = pdb.get_histogram();
        auto res = map.fit(pdb_h);
        fitter::FitReporter::report(res);
        fitter::FitReporter::save(res, path + "report.txt");

        plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
        plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

        auto scan = map.cutoff_scan(100, pdb_h);
        plots::PlotDataset::quick_plot(scan.as_dataset(), path + "scan.pdf");
    }

    //* MEASUREMENT FIT
    // Fit the measurements to the EM density map.
    if (true) {
        std::string path = "figures/em/measurement_fit/" + utility::stem(mfile) + "/" + utility::stem(mapfile) + "/";

        auto res = map.fit(mfile);
        fitter::FitReporter::report(res);
        fitter::FitReporter::save(res, path + "report.txt");

        plots::PlotIntensityFit::quick_plot(res, path + "intensity_fit.pdf");
        plots::PlotIntensityFitResiduals::quick_plot(res, path + "residuals.pdf");

        // auto scan = map.cutoff_scan(100, mfile);
        // plots::PlotDataset::quick_plot(scan, path + "scan.pdf");
    }

    //* GENERATE INTENSITY & PDB
    // Perform an intensity fit, and generate a PDB file for the optimal parameters.
    if (false) {
        // unsigned int c;
        // for(const em::Image& image : map.images()) {
        //     plots::PlotImage::quick_plot(image, "figures/em/images/" + utility::stem(mapfile) + "/" + std::to_string(c++) + ".png");
        // }
        std::string path = "figures/em/cut/" + utility::stem(mapfile) + "/";
        // plots::PlotIntensity::quick_plot(map.get_histogram(map.level(3)), path + "intensity3.pdf");
        // plots::PlotIntensity::quick_plot(map.get_histogram(map.level(2.5)), path + "intensity25.pdf");
        // plots::PlotIntensity::quick_plot(map.get_histogram(map.level(2)), path + "intensity2.pdf");
        // plots::PlotIntensity::quick_plot(map.get_histogram(map.level(1.5)), path + "intensity15.pdf");
        // plots::PlotIntensity::quick_plot(map.get_histogram(map.level(1)), path + "intensity1.pdf");

        auto d1 = map.get_histogram(map.from_level(3)).calc_debye_scattering_intensity();
        d1.normalize(1e7);
        plots::PlotIntensity plot(d1);
        auto p = map.get_protein(map.from_level(3));
        p->generate_new_hydration();
        p->get_histogram().calc_debye_scattering_intensity().save("hydrated.txt");

        auto d = p->get_histogram().calc_debye_scattering_intensity();
        d.normalize(1e7);
        plot.plot(d);
        plot.save("test.pdf");
        // plots::PlotIntensity::quick_plot(p->get_histogram(), "test.pdf");

        // map.get_histogram(map.level(3)).calc_debye_scattering_intensity().save("intensity3.txt");
        // map.get_histogram(map.level(2.5)).calc_debye_scattering_intensity().save("intensity25.txt");
        // map.get_histogram(map.level(2)).calc_debye_scattering_intensity().save("intensity2.txt");
        // map.get_histogram(map.level(1.5)).calc_debye_scattering_intensity().save("intensity15.txt");
        // map.get_histogram(map.level(1)).calc_debye_scattering_intensity().save("intensity1.txt");
        map.save(map.from_level(3), "data/output/" + utility::stem(mapfile) + ".pdb");
    }

    return 0;
}