#include "fitter/ExcludedVolumeFitter.h"
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
#include <hist/distance_calculator/HistogramManagerMTFFGridSurface.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

#include <array>

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb> <saxs>" << std::endl;
        return 1;
    }
    io::ExistingFile pdb(argv[1]);
    io::ExistingFile saxs(argv[2]);

    settings::general::output += "vary_grid_radii/" + pdb.stem() + "/";
    settings::axes::qmin = 1e-2;
    settings::axes::qmax = 1;
    settings::molecule::use_effective_charge = false;
    settings::grid::cell_width = 0.1;
    settings::grid::exv::surface_thickness = 1;

    data::Molecule molecule(pdb);
    molecule.generate_new_hydration();
    Dataset standard, profile_atomic;
    {
        auto hist = hist::HistogramManagerMTFFGridSurface(molecule).calculate_all();
        standard = hist->debye_transform().as_dataset();
        profile_atomic = hist->get_profile_aa().as_dataset();
    }
    settings::general::verbose = false;

    std::vector<Dataset> datasets, datasets_hydro_fitted, datasets_exv_fitted, profile_exv;
    std::vector<double> chi2_hydro, chi2_exv;
    for (double r = 1; r <= 3; r += 0.1) {
        std::cout << "radius: " << r << std::endl;
        settings::grid::min_exv_radius = r;
        auto hist = hist::HistogramManagerMTFFGridSurface(molecule).calculate_all();
        auto hist_cast = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get());
        datasets.push_back(hist->debye_transform().as_dataset());
        profile_exv.push_back(hist_cast->get_profile_xx().as_dataset());

        fitter::HydrationFitter fitter(saxs, std::move(hist));
        auto fit = fitter.fit();
        datasets_hydro_fitted.push_back(fit->info.fitted_intensity_interpolated);
        chi2_hydro.push_back(fit->fval/fit->dof);

        fitter::ExcludedVolumeFitter fitter_exv(saxs, hist::HistogramManagerMTFFGridSurface(molecule).calculate_all());
        auto fit_exv = fitter_exv.fit();
        datasets_exv_fitted.push_back(fit_exv->info.fitted_intensity_interpolated);
        chi2_exv.push_back(fit_exv->fval/fit_exv->dof);

        // molecule.get_grid()->save(settings::general::output + "grid_" + std::to_string(r) + ".pdb");
        molecule.clear_grid();
    }

    std::array<double, 3> rgb_start = {153, 0, 0}, rgb_end = {0, 0, 153};
    plots::PlotDataset plot, plot_hydro, plot_exv, plot_relative, plot_aa_xx;
    for (unsigned int i = 0; i < datasets.size(); i++) {
        double r = 1 + 0.1 * i;
        double r_norm = (r - 1) / 2;
        std::array<double, 3> rgb = {rgb_start[0] * (1 - r_norm) + rgb_end[0] * r_norm, rgb_start[1] * (1 - r_norm) + rgb_end[1] * r_norm, rgb_start[2] * (1 - r_norm) + rgb_end[2] * r_norm};
        std::stringstream ss;
        ss << "#" << std::hex << std::setfill('0') << std::setw(2) << static_cast<int>(rgb[0]) << std::setw(2) << static_cast<int>(rgb[1]) << std::setw(2) << static_cast<int>(rgb[2]);
        plot.plot(
            datasets[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, 
                {"logx", true}, {"logy", true}, {"color", ss.str()}, {"normalize", true}, {"title", "Unfitted profiles"}}
            )
        );
        plot_hydro.plot(
            datasets_hydro_fitted[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, 
                {"logx", true}, {"logy", true}, {"color", ss.str()}, {"normalize", true}, {"legend", std::to_string(chi2_hydro[i])},
                {"title", "Fitted hydration"}}
            )
        );
        plot_exv.plot(
            datasets_exv_fitted[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, 
                {"logx", true}, {"logy", true}, {"color", ss.str()}, {"normalize", true}, {"legend", std::to_string(chi2_exv[i])},
                {"title", "Fitted exv and hydration"}}
            )
        );

        SimpleDataset relative = standard;
        SimpleDataset aa_xx = profile_atomic;
        for (unsigned int j = 0; j < relative.size(); j++) {
            relative.y(j) = datasets[i].y(j) / standard.y(j);
            aa_xx.y(j) /= profile_exv[i].y(j);
        }
        relative.normalize();
        aa_xx.normalize();
        plot_relative.plot(
            relative, 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "relative deviation"}, {"yrange", std::vector{0.5, 1.5}}, {"logx", true}, {"color", ss.str()}, 
                {"normalize", true}, {"title", "Relative deviation from standard"}, {"legend", std::to_string(r)}}
            )
        );
        plot_aa_xx.plot(
            aa_xx, 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "aa / xx"}, {"yrange", std::vector{0.5, 1.5}}, {"logx", true}, {"color", ss.str()}, 
                {"normalize", true}, {"title", "Vacuum over exv"}, {"legend", std::to_string(r)}}
            )
        );
    }
    plot.save(settings::general::output + "analysis.png");
    plot_hydro.save(settings::general::output + "analysis_hydro_fitted.png");
    plot_exv.save(settings::general::output + "analysis_exv_fitted.png");
    plot_aa_xx.save(settings::general::output + "analysis_aa_xx.png");
    plot_relative.save(settings::general::output + "analysis_relative.png");
}