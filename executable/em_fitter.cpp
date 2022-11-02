#include <CLI/CLI.hpp>

#include <iostream>

#include <plots/all.h>
#include <em/ImageStack.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>

using std::string;

int main(int argc, char const *argv[]) {
    setting::protein::use_effective_charge = false;
    setting::em::hydrate = true;
    setting::fit::verbose = false;

    CLI::App app{"Fit an EM map to a SAXS measurement."};

    std::string mfile, mapfile, settings, output;
    app.add_option("input_map", mapfile, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input_exp", mfile, "Path to the SAXS measurement.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", output, "Path to save the generated figures at.");
    app.add_option("--qlow", setting::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qhigh", setting::axes::qmax, "Upper limit on used q values from measurement file.");
    app.add_option("--frequency", setting::em::sample_frequency, "Sampling frequency of the EM map.");
    CLI11_PARSE(app, argc, argv);

    // if a settings file was provided
    if (p_settings->count() != 0) {
        setting::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        setting::discover(std::filesystem::path(mfile).parent_path().string());
    }

    std::cout << setting::em::sample_frequency << std::endl;

    if (output.empty()) {
        output = "figures/";
    }
    setting::plot::path = output + "em_fitter/" + utility::stem(mfile) + "/" + utility::stem(mapfile) + "/";
    
    std::cout << "Performing EM fit with map " << mapfile << " and measurement " << mfile << std::endl;
    em::ImageStack map(mapfile); 

    // Fit the measurements to the EM density map.
    auto res = map.fit(mfile);
    FitReporter::report(res);
    FitReporter::save(setting::plot::path + "report.txt", res);

    res->figures.data.save(setting::plot::path + utility::stem(mfile) + ".dat");
    res->figures.intensity_interpolated.save(setting::plot::path + "fit.fit");
    plots::PlotIntensityFit::quick_plot(res, setting::plot::path + "intensity_fit.pdf");
    plots::PlotIntensityFitResiduals::quick_plot(res, setting::plot::path + "residuals.pdf");
    return 0;
}