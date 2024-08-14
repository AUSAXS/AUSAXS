#include "fitter/HydrationFitter.h"
#include "fitter/LinearFitter.h"
#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

#include <array>

int main(int argc, char const *argv[]) {
    io::ExistingFile pdb(argv[1]);
    io::ExistingFile saxs(argv[2]);

    //### GENERATE INTERNAL PLOT ###//
    settings::general::output += "vary_grid_radii/" + pdb.stem() + "/";
    settings::axes::qmin = 1e-2;
    settings::axes::qmax = 1;
    settings::molecule::use_effective_charge = false;
    settings::grid::cell_width = 0.1;

    data::Molecule molecule(pdb);
    molecule.generate_new_hydration();
    Dataset standard = hist::HistogramManagerMTFFGrid(molecule).calculate_all()->debye_transform().as_dataset();

    std::vector<Dataset> datasets;
    std::vector<Dataset> datasets_fitted;
    std::vector<double> chi2;
    for (double r = 1; r <= 3; r += 0.1) {
        settings::grid::min_exv_radius = r;
        auto hist = hist::HistogramManagerMTFFGrid(molecule).calculate_all();
        datasets.push_back(hist->debye_transform().as_dataset());
        fitter::HydrationFitter fitter(saxs, std::move(hist));
        auto fit = fitter.fit();
        datasets_fitted.push_back(fit->info.fitted_intensity_interpolated);
        chi2.push_back(fit->fval/fit->dof);
        // molecule.get_grid()->save(settings::general::output + "grid_" + std::to_string(r) + ".pdb");
        molecule.clear_grid();
    }

    std::array<double, 3> rgb_start = {153, 0, 0}, rgb_end = {0, 0, 153};
    plots::PlotDataset plot, plot_fitted, plot_relative;
    for (unsigned int i = 0; i < datasets.size(); i++) {
        double r = 1 + 0.1 * i;
        double r_norm = (r - 1) / 2;
        std::array<double, 3> rgb = {rgb_start[0] * (1 - r_norm) + rgb_end[0] * r_norm, rgb_start[1] * (1 - r_norm) + rgb_end[1] * r_norm, rgb_start[2] * (1 - r_norm) + rgb_end[2] * r_norm};
        std::stringstream ss;
        ss << "#" << std::hex << std::setfill('0') << std::setw(2) << static_cast<int>(rgb[0]) << std::setw(2) << static_cast<int>(rgb[1]) << std::setw(2) << static_cast<int>(rgb[2]);
        plot.plot(datasets[i], plots::PlotOptions({{"xlabel", "Distance"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, {"logx", true}, {"logy", true}, {"color", ss.str()}, {"normalize", true}}));
        plot_fitted.plot(datasets_fitted[i], plots::PlotOptions({{"xlabel", "Distance"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, {"logx", true}, {"logy", true}, {"color", ss.str()}, {"normalize", true}, {"legend", std::to_string(chi2[i])}}));

        SimpleDataset relative = standard;
        for (unsigned int j = 0; j < relative.size(); j++) {
            relative.y(j) = datasets[i].y(j) / standard.y(j);
        }
        relative.normalize();
        plot_relative.plot(relative, plots::PlotOptions({{"xlabel", "Distance"}, {"ylabel", "Intensity"}, {"yrange", std::vector{0.5, 1.5}}, {"logx", true}, {"color", ss.str()}, {"normalize", true}}));
    }
    plot.save(settings::general::output + "analysis.png");
    plot_fitted.save(settings::general::output + "analysis_fitted.png");
    plot_relative.save(settings::general::output + "analysis_relative.png");
}