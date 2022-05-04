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
    plots::PlotIntensity plot(data);

    // prepare fit colors
    gStyle->SetPalette(kSolar);
    auto cols = TColor::GetPalette();

    vector<Fit> fits;
    vector<string> paths;
    unsigned int color_step = (cols.GetSize()-1)/argc;
    Dataset fitted_vals("resolution", "cutoff");
    for (int i = 2; i < argc; i++) {
        std::cout << "Now fitting " << argv[i] << "..." << std::endl;
        string current_file = argv[i];
        em::ImageStack image(current_file);
        auto hist = protein.get_histogram();
        auto fit = image.fit(hist);

        // prepare writing the fit results
        fits.push_back(*fit);
        paths.push_back(current_file);

        // prepare dataset for concentration v. resolution plot
        fitted_vals.x.push_back(utility::extract_number<unsigned int>(current_file));
        fitted_vals.y.push_back(fit->params["cutoff"]);
        fitted_vals.yerr.push_back(fit->errors["cutoff"]);

        // actually plot this iteration
        plot.plot_intensity(fit, cols.At(color_step*(i-1)));

        // chi2 landscape plot
        Dataset contour = image.cutoff_scan({10, 0, 6}, hist);
        Dataset evaluated_points = image.fit(hist)->evaluated_points;
        evaluated_points.plot_options.set("markers", {{"color", kOrange+2}});

        plots::PlotDataset plot(contour);
        plot.plot(evaluated_points);
        plot.save("figures/stuff/landscapes/" + std::filesystem::path(current_file).stem().string() + ".pdf");
    }
    
    // write fit report to disk
    FitReporter::report(fits, paths);
    FitReporter::save("figures/fits/EMfit.txt", fits, paths);

    // generate concentration v. resolution plot
    fitted_vals.set_plot_options("errors");
    plots::PlotDataset d_plot(fitted_vals);
    d_plot.save("figures/stuff/fitted_vals.pdf");

    plot.save("figures/stuff/fits.pdf");
    return 0;
}