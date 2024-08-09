#include <CLI/CLI.hpp>

#include <data/record/Water.h>
#include <data/Molecule.h>
#include <fitter/HydrationFitter.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <io/ExistingFile.h>
#include <constants/Constants.h>
#include <form_factor/DisplacedVolumeTable.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <plots/All.h>
#include <settings/All.h>

#include <vector>
#include <string>
#include <iostream>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());
auto gaussian_noise = [] (double stddev) {
    std::normal_distribution<double> d(0, stddev);
    return d(gen);
};

auto generate_voronoi_exv_table = [] () {
    constants::displaced_volume::detail::DisplacedVolumeSet unc {
        .H   = 0,
        .C   = 0.9,
        .CH  = 1.4,
        .CH2 = 2.2,
        .CH3 = 3.5,
        .N   = 0.6,
        .NH  = 1.8,
        .NH2 = 3.8,
        .NH3 = 3.1,
        .O   = 3.0,
        .OH  = 4.0,
        .S   = 4.9,
        .SH  = 4.1
    };
    auto table = constants::displaced_volume::Voronoi_implicit_H;
    table.C   += gaussian_noise(unc.C);
    table.CH  += gaussian_noise(unc.CH);
    table.CH2 += gaussian_noise(unc.CH2);
    table.CH3 += gaussian_noise(unc.CH3);
    table.N   += gaussian_noise(unc.N);
    table.NH  += gaussian_noise(unc.NH);
    table.NH2 += gaussian_noise(unc.NH2);
    table.NH3 += gaussian_noise(unc.NH3);
    table.O   += gaussian_noise(unc.O);
    table.OH  += gaussian_noise(unc.OH);
    table.S   += gaussian_noise(unc.S);
    table.SH  += gaussian_noise(unc.SH);

    form_factor::storage::detail::set_custom_displaced_volume_table(table);
};

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    io::ExistingFile pdb, mfile;

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    CLI11_PARSE(app, argc, argv);

    settings::general::output = "output/vary_exv_tables/";
    settings::molecule::throw_on_unknown_atom = false;
    settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Custom;
    settings::molecule::use_effective_charge = false;
    settings::hist::fit_excluded_volume = true;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;

    unsigned int max = 1000;
    data::Molecule molecule(pdb);
    SimpleDataset data(mfile);
    molecule.add_implicit_hydrogens();
    molecule.generate_new_hydration();
    auto hist = molecule.get_histogram();
    std::vector<double> chi2s;

    auto make_fit = [&] () {
        fitter::ExcludedVolumeFitter fitter(
            data, 
            std::make_unique<hist::CompositeDistanceHistogramFFExplicit>(*static_cast<hist::CompositeDistanceHistogramFFExplicit*>(hist.get()))
        );
        auto fit = fitter.fit();
        return fit->fval/fit->dof;
    };

    for (unsigned int i = 0; i < max; i++) {
        std::cout << 100*double(i)/max << "%\r" << std::flush;
        generate_voronoi_exv_table();
        chi2s.push_back(make_fit());
    }
    std::cout << "min is " << *std::min_element(chi2s.begin(), chi2s.end()) << std::endl;

    Axis axis(0, 1000, 10000);
    auto bins = axis.as_vector();
    std::vector<double> chi2_hist(bins.size(), 0);
    for (auto chi2 : chi2s) {
        chi2_hist[std::min<int>(axis.get_bin(chi2), bins.size()-1)]++;
    }
    std::transform(chi2_hist.begin(), chi2_hist.end(), chi2_hist.begin(), [&] (double x) {return x/max;});
    double count = 0, end = chi2_hist.size();
    while (count < 0.01 && 0 < end) {
        count += chi2_hist[--end];
    }
    if (100 < end) {
        // merge histogram bins to get under 100
        unsigned int merge = end/100;
        std::vector<double> new_chi2_hist;
        for (unsigned int i = 0; i < chi2_hist.size()-merge; i += merge) {
            new_chi2_hist.push_back(std::accumulate(chi2_hist.begin()+i, chi2_hist.begin()+i+merge, 0.0));
        }
        chi2_hist = new_chi2_hist;
        bins = Axis(0, 1000, chi2_hist.size()).as_vector();
    }

    settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Voronoi_implicit_H;
    auto chi2_mean = make_fit();

    settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Traube;
    auto chi2_traube = make_fit();

    end = std::max<int>(end, axis.get_bin(chi2_traube)+10);
    end = std::max<int>(end, axis.get_bin(chi2_mean)+10);

    plots::PlotDataset plot(SimpleDataset{bins, chi2_hist}, plots::PlotOptions({{"xlabel", "$\\chi^2_r$"}, {"ylabel", "density"}, {"xlim", Limit{0, axis.get_bin_value(end)}}}));
    plot.vline(chi2_mean, plots::PlotOptions({{"color", "k"}, {"legend", "mean"}, {"ls", style::line::dashed}}));
    plot.vline(chi2_traube, plots::PlotOptions({{"color", "r"}, {"legend", "Traube"}, {"ls", style::line::dashed}}));
    plot.save(settings::general::output + mfile.stem() + ".png");
    return 0;
}