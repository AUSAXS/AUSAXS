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

    settings::general::output += "vary_grid_exv_density/" + pdb.stem() + "/";
    settings::molecule::use_effective_charge = false;
    settings::axes::qmax = 1;

    data::Molecule molecule(pdb);
    molecule.generate_new_hydration();

    auto hist = hist::HistogramManagerMTFFGridSurface(molecule).calculate_all();
    SimpleDataset profile_vacuum = hist->get_profile_aa().as_dataset();
    std::vector<SimpleDataset> datasets, profile_exv;

    Axis scale_axis = {0.9, 1.1, 20};
    for (double scale = scale_axis.min; scale < scale_axis.max; scale += scale_axis.step()) {
        // form_factor::ExvFormFactor ffx(std::pow(settings::grid::exv::width, 3)*scale);
        
        auto V = std::pow(settings::grid::exv::width, 3);
        form_factor::ExvFormFactor ffx(V);
        ffx.q0 *= scale;
        hist::CompositeDistanceHistogramFFGrid::regenerate_ff_table(std::move(ffx));

        auto hist_cast = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get());
        hist_cast->apply_excluded_volume_scaling_factor(1+scale*1e-6); // trigger internal recalculation
        datasets.push_back(hist.get()->debye_transform().as_dataset());
        profile_exv.push_back(hist_cast->get_profile_xx().as_dataset());
    }

    plots::PlotDataset plot, plot_aa_xx, plot_xx;
    style::color_map::Rainbow color_map(datasets.size());
    for (unsigned int i = 0; i < datasets.size(); i++) {
        auto c = color_map.next();
        plot.plot(
            datasets[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"xrange", std::vector{1e-2, 1.}}, 
                {"logx", true}, {"logy", true}, {"color", c}, {"normalize", true}, {"title", "Unfitted profiles"},
                {"legend", "$\\sigma = $" + std::to_string(scale_axis.get_bin_value(i))}}
            )
        );

        plot_xx.plot(
            profile_exv[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"xrange", std::vector{1e-2, 1.}}, 
                {"logx", true}, {"logy", true}, {"color", c}, {"title", "Unfitted profiles"},
                {"legend", "$\\sigma = $" + std::to_string(scale_axis.get_bin_value(i))}}
            )
        );

        SimpleDataset aa_xx = profile_exv[i];
        for (unsigned int j = 0; j < profile_exv[i].size(); ++j) {
            aa_xx.y(j) /= profile_vacuum.y(j);
        }
        plot_aa_xx.plot(
            aa_xx, 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "aa / xx"}, {"logx", true}, {"color", c}, {"title", "Vacuum over exv"}, 
                {"legend", "$\\sigma = $" + std::to_string(scale_axis.get_bin_value(i))}, {"xrange", std::vector{1e-2, 1.}}}
            )
        );
    }
    plot.save(settings::general::output + "profiles.png");
    plot_aa_xx.save(settings::general::output + "aa_xx.png");
    plot_xx.save(settings::general::output + "xx.png");
}