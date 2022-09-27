#include <CLI/CLI.hpp>

#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    setting::em::hydrate = true;

    CLI::App app{"Fit an EM map to a SAXS measurement."};

    std::string mfile, mapfile, settings, output;
    app.add_option("input_map", mapfile, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input_exp", mfile, "Path to the SAXS measurement.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", output, "Path to save the generated figures at.");
    app.add_option("--qlow", setting::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qhigh", setting::axes::qmax, "Upper limit on used q values from measurement file.");
    CLI11_PARSE(app, argc, argv);

    if (p_settings->count() != 0) {
        setting::read(settings);
        CLI11_PARSE(app, argc, argv);
    } else {
        std::string path = std::filesystem::path(mfile).parent_path().string();
        if (std::filesystem::exists(path + "/settings.txt")) {
            std::cout << "Using discovered settings file at " << path << "/settings.txt" << std::endl;
            setting::read(path + "/settings.txt");
        }
    }

    if (output.empty()) {
        output = "figures/em_fitter/" + utility::stem(mfile) + "/" + utility::stem(mapfile) + "/";
    }
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