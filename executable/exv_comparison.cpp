#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <settings/All.h>
#include <plots/All.h>

int main(int argc, char const *argv[]) {
    std::string s_pdb, s_mfile;
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", s_pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", s_mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/intensity_fitter/")->group("General options");

    data::Molecule molecule(s_pdb);
    auto mtffg   = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFGrid<false>(molecule).calculate_all().get())->get_profile_xx();
    auto mtffavg = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFAvg<false>(molecule).calculate_all().get())->get_profile_xx();
    auto mtffexp = static_cast<hist::ICompositeDistanceHistogramExv*>(hist::HistogramManagerMTFFExplicit<false>(molecule).calculate_all().get())->get_profile_xx();

    plots::PlotIntensity();
}