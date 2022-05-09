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
    setting::em::sample_frequency = 1;
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

        // chi2 contour plot
        Dataset contour = image.cutoff_scan({100, 0, 6}, hist);
        Dataset evaluated_points = fit->evaluated_points;
        evaluated_points.plot_options.set("markers", {{"color", kOrange+2}});

        plots::PlotDataset plot_c(contour);
        plot_c.plot(evaluated_points);
        plot_c.save("figures/stuff/landscapes/" + std::filesystem::path(current_file).stem().string() + ".png");

        // generate residual v. resolution plot
        Dataset& residuals = fit->residuals;
        residuals.add_plot_options("errors", {{"logx", true}, {"xlabel", "q"}, {"ylabel", "residual"}});
        plots::PlotDataset::quick_plot(fitted_vals, "figures/stuff/residuals/" + std::filesystem::path(current_file).stem().string() + ".pdf");
    }
    
    // write fit report to disk
    FitReporter::report(fits, paths);
    FitReporter::save("figures/fits/EMfit.txt", fits, paths);

    // generate concentration v. resolution plot
    fitted_vals.add_plot_options("errors");
    plots::PlotDataset::quick_plot(fitted_vals, "figures/stuff/cutoff.pdf");

    // generate chi2 v. resolution plot
    vector<double> chi2(fits.size());
    std::transform(fits.begin(), fits.end(), chi2.begin(), [] (const Fit& fit) {return fit.chi2;});
    fitted_vals.y = chi2;
    fitted_vals.add_plot_options({{"ylabel", "chi2"}});
    plots::PlotDataset::quick_plot(fitted_vals, "figures/stuff/chi2.pdf");

    // generate intensity comparison plots
    Multiset intensities;
    std::transform(fits.begin(), fits.end(), std::back_inserter(intensities.data), [] (const Fit& fit) {return fit.figures[1];});
    plots::PlotResolutionComparison plot_r(intensities);
    plot_r.save("figures/stuff/fits.pdf");
    return 0;
}