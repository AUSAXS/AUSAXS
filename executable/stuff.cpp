#include <iostream>
#include <filesystem>

#include <em/ImageStack.h>
#include <plots/all.h>
#include <fitter/FitReporter.h>
#include <utility/Exceptions.h>
#include <utility/Settings.h>
#include <utility/Utility.h>
#include <math/Statistics.h>

using std::string;

int main(int argc, char const *argv[]) {
    unsigned int loops = 2; // how many times each map is fitted to get an average chi2

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
    Dataset voxel_sizes("resolution", "voxel size [Angstrom]");
    for (int i = 2; i < argc; i++) {
        std::cout << "Now fitting " << argv[i] << "..." << std::endl;
        string current_file = argv[i];
        unsigned int resolution = utility::extract_number<unsigned int>(current_file);
        em::ImageStack image(current_file, resolution);
        auto hist = protein.get_histogram();


        //#######################################//
        //#    Repeat fits & combine results    #//
        //#######################################//
        vector<Fit> cur_fits;
        for (unsigned int i = 0; i < loops; i++) {
            cur_fits.push_back(*image.fit(hist));
        }

        std::vector<double> holder(loops);
        Fit fit = cur_fits[0];

        // use mean chi2
        std::transform(cur_fits.begin(), cur_fits.end(), holder.begin(), [] (const Fit& fit) {return fit.chi2;});
        fit.chi2 = stats::mean(holder);

        // use total fevals
        fit.calls = std::accumulate(cur_fits.begin(), cur_fits.end(), 0, [] (unsigned int sum, const Fit& fit) {return sum + fit.calls;});

        // use mean cutoff
        stats::MeasurementSeries ms(loops);
        std::transform(cur_fits.begin(), cur_fits.end(), ms.begin(), [] (const Fit& fit) {return fit.get_parameter("cutoff");});

        auto& p_cutoff = fit.get_parameter("cutoff");
        p_cutoff.value = ms.weighted_mean();
        p_cutoff.set_error(ms.error_of_mean());

        fits.push_back(fit);
        paths.push_back(current_file);

        //#######################################//
        //#           Prepare plots             #//
        //#######################################//

        // prepare dataset for cutoff v. resolution plot
        auto cutoff = fit.get_parameter("cutoff");
        fitted_vals.push_back({double(resolution), cutoff.value, cutoff.mean_error()});

        // prepare voxel size v. resolution plot
        auto header = image.get_header();
        if (!utility::equal(header->cella_x, header->cella_y, header->cella_z) || 
            !utility::equal(header->nx, header->ny, header->nz)) {
                throw except::size_error("Error in stuff: Header dimensions are not equal!");
        }
        voxel_sizes.push_back({double(resolution), header->cella_x/header->nx});

        // chi2 contour plot
        // Dataset evaluated_points = fit.evaluated_points;
        // evaluated_points.plot_options.set("markers", {{"color", kOrange+2}});

        // plots::PlotDataset plot_c(contour);
        // plot_c.plot(evaluated_points);
        // plot_c.save("figures/stuff/landscapes/" + std::filesystem::path(current_file).stem().string() + ".png");

        // generate residual v. resolution plot
        Dataset& residuals = fit.residuals;
        residuals.add_plot_options("errors", {{"logx", true}, {"xlabel", "q"}, {"ylabel", "residual"}});
        plots::PlotDataset::quick_plot(residuals, "figures/stuff/residuals/" + std::filesystem::path(current_file).stem().string() + ".png");
    }

    // write fit report to disk
    FitReporter::report(fits, paths);
    FitReporter::save("figures/fits/EMfit.txt", fits, paths);

    // generate cutoff v. resolution plot
    fitted_vals.add_plot_options("markers");
    plots::PlotDataset::quick_plot(fitted_vals, "figures/stuff/cutoff.pdf");

    // generate chi2 v. resolution plot
    vector<double> chi2(fits.size());
    std::transform(fits.begin(), fits.end(), chi2.begin(), [] (const Fit& fit) {return fit.chi2;});
    fitted_vals.y = chi2;
    fitted_vals.add_plot_options({{"ylabel", "chi2"}});
    plots::PlotDataset::quick_plot(fitted_vals, "figures/stuff/chi2.pdf");

    // generate voxel size v. resolution plot
    voxel_sizes.add_plot_options("marker");
    plots::PlotDataset::quick_plot(voxel_sizes, "figures/stuff/voxel_sizes.pdf");

    // generate intensity comparison plots
    Multiset intensities;
    std::transform(fits.begin(), fits.end(), std::back_inserter(intensities.data), [] (const Fit& fit) {return fit.figures.intensity;});
    intensities.save("temp/multiset");
    intensities.ylimits(1e-4, inf);
    plots::PlotResolutionComparison plot_r(intensities);
    plot_r.save("figures/stuff/fits.pdf");
    return 0;
}