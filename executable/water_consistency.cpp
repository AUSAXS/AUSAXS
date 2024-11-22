#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <fitter/SmartFitter.h>
#include <dataset/SimpleDataset.h>
#include <utility/Utility.h>
#include <settings/All.h>
#include <plots/All.h>

#include <fstream>
#include <iomanip>

using namespace ausaxs;

int main(int, char const*[]) {
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;
    settings::molecule::use_effective_charge = false;
    settings::fit::fit_exv_debye_waller = true;

    SimpleDataset saxs("data/hub_data/193l_lyz/193l_lyz.dat");
    data::Molecule protein("data/hub_data/193l_lyz/193l_lyz.pdb");
    protein.add_implicit_hydrogens();

    unsigned int N = 50; 
    std::vector<SimpleDataset> ww;
    for (unsigned int i = 0; i < N; ++i) {
        std::cout << "Progress: " << double(i)/N*100 << "%\r" << std::flush;
        protein.generate_new_hydration();

        fitter::SmartFitter fitter(saxs, protein.get_histogram());
        auto result = fitter.fit();
        ww.push_back(static_cast<hist::ICompositeDistanceHistogram*>(fitter.get_model())->get_profile_ww());
    }

    plots::PlotDataset plot(ww[0], plots::PlotOptions({{"logx", true}, {"logy", true}, {"xlabel", "q"}, {"ylabel", "I(q)"}}));
    for (unsigned int i = 1; i < N; ++i) {
        plot.plot(ww[i], plots::PlotOptions({{"color", style::color::blue}}));
    }
    plot.save("output/water_consistency/ww_profiles.png");

    return 0;
}