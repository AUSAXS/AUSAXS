#include <iostream>
#include <filesystem>

#include <em/ImageStack.h>
#include <plots/all.h>
#include <fitter/FitReporter.h>
#include <utility/Exceptions.h>
#include <utility/Settings.h>
#include <utility/Utility.h>

using std::string;

int main(int argc, char const *argv[]) {
    string pdb_file = argv[1];

    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    setting::fit::q_high = 0.4;

    // plot pdb file
    Protein protein(pdb_file);
    SAXSDataset data = protein.get_histogram().calc_debye_scattering_intensity();
    data.reduce(setting::fit::N, true);
    data.limit(Limit(setting::fit::q_low, setting::fit::q_high));
    data.simulate_errors();

    vector<Fit> fits;
    vector<string> paths;
    Dataset fitted_vals("resolution", "cutoff");
    Multiset intensities;
    for (int i = 2; i < argc; i++) {
        std::cout << "Now fitting " << argv[i] << "..." << std::endl;
        string current_file = argv[i];
        unsigned int resolution = utility::extract_number<unsigned int>(current_file);
        em::ImageStack image(current_file, resolution);
        auto hist = protein.get_histogram();
        auto fit = image.fit(hist);

        // prepare writing the fit results
        fits.push_back(*fit);
        paths.push_back(current_file);

        // prepare dataset for cutoff v. resolution plot
        fitted_vals.x.push_back(resolution);
        fitted_vals.y.push_back(fit->params["cutoff"]);
        fitted_vals.yerr.push_back(fit->errors["cutoff"]);

        // actually plot this iteration
        intensities.push_back(fit->figures[1]);

        // chi2 contour plot
        // Dataset contour = image.cutoff_scan({100, 0, 6}, hist);
        // Dataset evaluated_points = fit->evaluated_points;
        // evaluated_points.plot_options.set("markers", {{"color", kOrange+2}});

        // plots::PlotDataset plot_c(contour);
        // plot_c.plot(evaluated_points);
        // plot_c.save("figures/stuff/landscapes/" + std::filesystem::path(current_file).stem().string() + ".pdf");
    }
    
    // write fit report to disk
    FitReporter::report(fits, paths);
    FitReporter::save("figures/fits/EMfit.txt", fits, paths);

    // generate concentration v. resolution plot
    fitted_vals.add_plot_options("errors");
    plots::PlotDataset d_plot(fitted_vals);
    d_plot.save("figures/stuff/fitted_vals.pdf");

    plots::PlotResolutionComparison plot_r(intensities);
    plot_r.save("figures/stuff/fits.pdf");
    return 0;
}