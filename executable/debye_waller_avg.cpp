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
#include <math/Statistics.h>

#include <random>

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb> <saxs>" << std::endl;
        return 1;
    }
    io::ExistingFile pdb(argv[1]);
    io::ExistingFile saxs(argv[2]);

    settings::general::output += "debye_waller_avg/" + pdb.stem() + "/";
    settings::molecule::use_effective_charge = false;
    settings::axes::qmax = 1;

    data::Molecule molecule(pdb);
    auto atoms = molecule.get_atoms();
    auto hist = hist::HistogramManagerMT<true>(molecule).calculate_all();
    SimpleDataset profile_vacuum = hist->get_profile_aa().as_dataset();

    auto stds = Axis(0.1, 1, 9).as_vector();
    unsigned int loops = 10;
    for (double std : stds) {
        auto add_noise = [std] (data::Molecule& mol) {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            static std::normal_distribution<> gauss(0, std);

            for (auto& body : mol.get_bodies()) {
                for (auto& atom : body.get_atoms()) {
                    atom.get_coordinates() += Vector3<double>(gauss(gen), gauss(gen), gauss(gen));
                }
            }
        };

        std::vector<Dataset> vacuums;
        for (unsigned int i = 0; i < loops; ++i) {
            data::Molecule molecule_copy(atoms);
            add_noise(molecule_copy);
            vacuums.push_back(hist::HistogramManagerMT<true>(molecule_copy).calculate_all()->get_profile_aa().as_dataset());
        }

        SimpleDataset average;
        for (unsigned int i = 0; i < profile_vacuum.size(); ++i) {
            std::vector<double> intensities;
            for (auto& vacuum : vacuums) {
                intensities.push_back(vacuum.y(i));
            }
            average.push_back(profile_vacuum.x(i), stats::mean(intensities), stats::std(intensities));
        }

        //################//
        //### Plotting ###//
        //################//
        plots::PlotDataset plot_vacuums;
        style::color_map::Rainbow color_map(vacuums.size());
        for (unsigned int i = 0; i < vacuums.size(); i++) {
            auto c = color_map.next();
            plot_vacuums.plot(
                vacuums[i], 
                plots::PlotOptions({
                    {"xlabel", "q"}, {"ylabel", "Intensity"}, {"xrange", std::vector{1e-2, 1.}}, 
                    {"logx", true}, {"logy", true}, {"color", c}, {"title", "profiles"}
                })
            );
        }
        plot_vacuums.save(settings::general::output + std::to_string(std) + "/vacuum_profiles.png");

        plots::PlotDataset plot_avg;
        {
            plot_avg.plot(
                average, 
                plots::PlotOptions({
                    {"xlabel", "q"}, {"ylabel", "Intensity"}, {"xrange", std::vector{1e-2, 1.}}, 
                    {"logx", true}, {"logy", true}, {"color", "black"}, {"title", "average"}
                })
            );

            auto max = average, min = average;
            for (unsigned int i = 0; i < average.size(); i++) {
                max.y(i) += average.yerr(i);
                min.y(i) -= average.yerr(i);
            }
            plot_avg.plot(max, plots::PlotOptions({{"color", "black"}, {"linestyle", "dashed"}}));
            plot_avg.plot(min, plots::PlotOptions({{"color", "black"}, {"linestyle", "dashed"}}));
        }

        plot_avg.save(settings::general::output + std::to_string(std) + "/average_profile.png");
    }
}