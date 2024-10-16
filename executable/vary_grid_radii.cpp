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
#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/HydrationFitter.h>
#include <fitter/LinearFitter.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb> <saxs>" << std::endl;
        return 1;
    }
    io::ExistingFile pdb(argv[1]);
    io::ExistingFile saxs(argv[2]);

    settings::general::output += "vary_grid_radii/" + pdb.stem() + "/";
    settings::molecule::use_effective_charge = false;
    settings::grid::cell_width = 0.5;
    settings::grid::exv::width = settings::grid::cell_width;
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

    std::vector<Dataset> datasets, datasets_hydro_fitted, profile_exv;
    std::vector<double> chi2_hydro, volumes, dummy_atoms;
    for (double r = 1.5; r <= 3; r += 0.1) {
        molecule.clear_grid();
        std::cout << "radius: " << r << std::endl;
        settings::grid::min_exv_radius = r;

        volumes.push_back(molecule.get_volume_grid());
        dummy_atoms.push_back(molecule.get_grid()->generate_excluded_volume(false).interior.size());

        auto hist = hist::HistogramManagerMTFFGridSurface(molecule).calculate_all();
        auto hist_cast = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get());
        datasets.push_back(hist->debye_transform().as_dataset());
        profile_exv.push_back(hist_cast->get_profile_xx().as_dataset());

        fitter::HydrationFitter fitter(saxs, std::move(hist));
        auto fit = fitter.fit();
        datasets_hydro_fitted.push_back(fit->info.fitted_intensity_interpolated);
        chi2_hydro.push_back(fit->fval/fit->dof);

        // molecule.get_grid()->save(settings::general::output + "grid_" + std::to_string(r) + ".pdb");
    }

    style::color_map::Rainbow cmap(datasets.size());
    plots::PlotDataset plot, plot_hydro, plot_relative, plot_aa_xx;
    std::vector<double> indices(datasets.size());
    for (unsigned int i = 0; i < datasets.size(); i++) {
        indices[i] = i;
        double r = 1.5 + 0.1 * i;
        auto c = cmap.next();
        plot.plot(
            datasets[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, 
                {"logx", true}, {"logy", true}, {"color", c}, {"normalize", true}, {"title", "Unfitted profiles"}}
            )
        );
        plot_hydro.plot(
            datasets_hydro_fitted[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, 
                {"logx", true}, {"logy", true}, {"color", c}, {"normalize", true}, {"legend", std::to_string(chi2_hydro[i])},
                {"title", "Fitted hydration"}}
            )
        );

        SimpleDataset relative = standard;
        SimpleDataset aa_xx = profile_exv[i];
        for (unsigned int j = 0; j < relative.size(); j++) {
            relative.y(j) = datasets[i].y(j) / standard.y(j);
            aa_xx.y(j) /= profile_atomic.y(j);
        }
        relative.normalize();
        aa_xx.normalize();
        plot_relative.plot(
            relative, 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "relative deviation"}, {"yrange", std::vector{0.5, 1.5}}, {"logx", true}, {"color", c}, 
                {"normalize", true}, {"title", "Relative deviation from standard"}, {"legend", std::to_string(r)}}
            )
        );
        plot_aa_xx.plot(
            aa_xx, 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "xx / aa"}, {"yrange", std::vector{0.5, 1.5}}, {"logx", true}, {"color", c}, 
                {"normalize", true}, {"title", "Vacuum over exv"}, {"legend", std::to_string(r)}}
            )
        );
    }
    plot.save(settings::general::output + "analysis.png");
    plot_hydro.save(settings::general::output + "analysis_hydro_fitted.png");
    plot_aa_xx.save(settings::general::output + "analysis_aa_xx.png");
    plot_relative.save(settings::general::output + "analysis_relative.png");
    plots::PlotDataset::quick_plot(
        SimpleDataset{indices, volumes}, 
        plots::PlotOptions({{"xlabel", "iteration"}, {"ylabel", "volume [Å$^-1$]"}, {"title", "Grid volume"}, {"color", "k"}}), 
        settings::general::output + "volumes.png"
    );
    plots::PlotDataset::quick_plot(
        SimpleDataset{indices, dummy_atoms}, 
        plots::PlotOptions({{"xlabel", "iteration"}, {"ylabel", "number of dummy exv atoms"}, {"title", "Exv dummy atoms"}, {"color", "k"}}), 
        settings::general::output + "exv_atoms.png"
    );
}