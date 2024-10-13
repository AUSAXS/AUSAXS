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

    settings::general::output += "vary_grid_ff/" + pdb.stem() + "/";
    settings::axes::qmin = 1e-2;
    settings::axes::qmax = 1;
    settings::molecule::use_effective_charge = false;

    data::Molecule molecule(pdb);
    molecule.generate_new_hydration();

    auto hist = hist::HistogramManagerMTFFGridSurface(molecule).calculate_all();
    Dataset profile_vacuum = hist->get_profile_aa().as_dataset();
    std::vector<Dataset> datasets, profile_exv;

    Axis scale_axis = {0.5, 1.5, 20};
    for (double scale = scale_axis.min; scale < scale_axis.max; scale += scale_axis.step()) {
        form_factor::ExvFormFactor ffx(std::pow(settings::grid::exv::width, 3)*scale);
        hist::CompositeDistanceHistogramFFGrid::regenerate_ff_table(std::move(ffx));

        // auto hist_cast = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get());
        // hist_cast->apply_excluded_volume_scaling_factor(1+scale*1e-6); // trigger internal recalculation
        // datasets.push_back(hist.get()->debye_transform().as_dataset());
        // profile_exv.push_back(hist_cast->get_profile_xx().as_dataset());

        auto h = hist::HistogramManagerMTFFGridSurface(molecule).calculate_all();
        auto h_cast = static_cast<hist::ICompositeDistanceHistogramExv*>(h.get());
        datasets.push_back(h->debye_transform().as_dataset());
        profile_exv.push_back(h_cast->get_profile_xx().as_dataset());
    }

    std::array<double, 3> rgb_start = {153, 0, 0}, rgb_end = {0, 0, 153};
    plots::PlotDataset plot, plot_aa_xx;
    for (unsigned int i = 0; i < datasets.size(); i++) {
        std::stringstream ss;
        {
            double r = 1 + 0.1 * i;
            double r_norm = (r - 1) / 2;
            std::array<double, 3> rgb = {rgb_start[0] * (1 - r_norm) + rgb_end[0] * r_norm, rgb_start[1] * (1 - r_norm) + rgb_end[1] * r_norm, rgb_start[2] * (1 - r_norm) + rgb_end[2] * r_norm};
            ss << "#" << std::hex << std::setfill('0') << std::setw(2) << static_cast<int>(rgb[0]) << std::setw(2) << static_cast<int>(rgb[1]) << std::setw(2) << static_cast<int>(rgb[2]);
        }
        plot.plot(
            datasets[i], 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "Intensity"}, {"yrange", std::vector{1e-4, 1e2}}, {"normalize", true}, 
                {"logx", true}, {"logy", true}, {"color", ss.str()}, {"normalize", true}, {"title", "Unfitted profiles"},
                {"legend", std::to_string(scale_axis.get_bin_value(i))}}
            )
        );

        SimpleDataset aa_xx = profile_vacuum;
        for (unsigned int j = 0; j < profile_exv.size(); ++j) {
            aa_xx.y(j) /= profile_exv[i].y(j);
        }
        aa_xx.normalize();
        plot_aa_xx.plot(
            aa_xx, 
            plots::PlotOptions(
                {{"xlabel", "q"}, {"ylabel", "aa / xx"}, {"yrange", std::vector{0.5, 1.5}}, {"logx", true}, {"color", ss.str()}, 
                {"normalize", true}, {"title", "Vacuum over exv"}, {"legend", std::to_string(scale_axis.get_bin_value(i))}}
            )
        );
    }
    plot.save(settings::general::output + "profiles.png");
    plot_aa_xx.save(settings::general::output + "aa_xx.png");
}