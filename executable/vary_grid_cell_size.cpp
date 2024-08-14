#include <CLI/CLI.hpp>

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <data/Molecule.h>
#include <settings/All.h>
#include <plots/All.h>

int main(int argc, char const *argv[]) {
    settings::general::output = "output/vary_grid_cell_size/";
    settings::molecule::use_effective_charge = false;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;
    settings::axes::qmax = 1;
    settings::grid::cell_width = 0.5;

    io::ExistingFile pdb(argv[1]);
    io::ExistingFile saxs(argv[2]);

    data::Molecule molecule(pdb);
    molecule.generate_new_hydration();

    std::vector<SimpleDataset> profiles;
    std::vector<double> exv_radii;
    Limit lims{0.25, 3.5};
    for (double exv_radius = lims.min; exv_radius <= lims.max; exv_radius += 0.25) {
        settings::grid::exv::radius = exv_radius;
        hist::CompositeDistanceHistogramFFGrid::regenerate_table();
        auto hist = molecule.get_histogram()->debye_transform();
        hist.normalize_max();
        profiles.push_back(std::move(hist));
        exv_radii.push_back(exv_radius);
    }

    std::vector<std::string> colors;
    {
        std::array<double, 3> rgb_start = {153, 0, 0}, rgb_end = {0, 0, 153};
        for (unsigned int i = 0; i < profiles.size(); i++) {
            double r = exv_radii[i];
            double r_norm = (r - lims.min) / lims.span();
            std::array<double, 3> rgb = {rgb_start[0] * (1 - r_norm) + rgb_end[0] * r_norm, rgb_start[1] * (1 - r_norm) + rgb_end[1] * r_norm, rgb_start[2] * (1 - r_norm) + rgb_end[2] * r_norm};
            std::stringstream ss;
            ss << "#" << std::hex << std::setfill('0') << std::setw(2) << static_cast<int>(rgb[0]) << std::setw(2) << static_cast<int>(rgb[1]) << std::setw(2) << static_cast<int>(rgb[2]);
            colors.push_back(ss.str());
        }
    }

    {   // absolute profiles
        plots::PlotDataset plot;
        for (unsigned int i = 0; i < profiles.size(); i++) {
            plot.plot(
                profiles[i], 
                plots::PlotOptions({
                    {"xlabel", "Distance"}, 
                    {"ylabel", "Intensity"}, 
                    {"yrange", std::vector{1e-3, 1.1}}, 
                    {"normalize", true}, 
                    {"logx", true}, 
                    {"logy", true}, 
                    {"normalize", true}, 
                    {"color", colors[i]},
                    {"legend", "$r_{exv}$ = " + std::to_string(exv_radii[i])}
                })
            );
        }
        plot.save(settings::general::output + "intensities.png");
    }

    {   // relative difference
        std::vector<SimpleDataset> diff = profiles;
        for (unsigned int i = 0; i < diff.size(); i++) {
            for (unsigned int j = 0; j < diff[i].size(); j++) {
                diff[i].y(j) = diff[i].y(j) / profiles[0].y(j);
            }
        }

        plots::PlotDataset plot;
        for (unsigned int i = 1; i < diff.size(); i++) {
            plot.plot(
                diff[i], 
                plots::PlotOptions({
                    {"xlabel", "Distance"}, 
                    {"ylabel", "Intensity"}, 
                    {"yrange", std::vector{0, 5}}, 
                    {"logx", true}, 
                    {"color", colors[i]},
                    {"legend", "$r_{exv}$ = " + std::to_string(exv_radii[i])}
                })
            );
        }
        plot.save(settings::general::output + "difference.png");
    }
}