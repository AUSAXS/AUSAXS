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
    static constants::displaced_volume::detail::DisplacedVolumeSet unc {
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
    return table;
};

auto generate_mf_exv_table = [] () {
    static constants::displaced_volume::detail::DisplacedVolumeSet unc {
        .H   = 0,
        .C   = 3.9,
        .CH  = 3.3,
        .CH2 = 1.2,
        .CH3 = 2.9,
        .N   = 0.4,
        .NH  = 1.7,
        .NH2 = 2.5,
        .NH3 = 2.6,
        .O   = 1.4,
        .OH  = 2.2,
        .S   = 6.6,
        .SH  = 9.0  
    };
    auto table = constants::displaced_volume::MinimumFluctuation_implicit_H;
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
    return table;
};

auto as_string = [] (const std::vector<double>& chi2, const std::vector<constants::displaced_volume::detail::DisplacedVolumeSet>& tables) {
    static auto table_to_string = [] (const constants::displaced_volume::detail::DisplacedVolumeSet& table) {
        std::stringstream ss;
        ss  << std::setprecision(3) << std::fixed
            << table.H << "\t" 
            << table.C << "\t" << table.CH << "\t" << table.CH2 << "\t" << table.CH3 << "\t" 
            << table.N << "\t" << table.NH << "\t" << table.NH2 << "\t" << table.NH3 << "\t" 
            << table.O << "\t" << table.OH << "\t" 
            << table.S << "\t" << table.SH;
        return ss.str();
    };

    static auto table_header_string = [] () {
        return "chi2\tH\tC\tCH\tCH2\tCH3\tN\tNH\tNH2\tNH3\tO\tOH\tS\tSH";
    };

    std::stringstream ss;
    ss << table_header_string() << std::endl;
    for (unsigned int i = 0; i < chi2.size(); i++) {
        ss << chi2[i] << "\t" << table_to_string(tables[i]) << std::endl;
    }
    return ss.str();
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

    unsigned int repeats = 1000;
    data::Molecule molecule(pdb);
    SimpleDataset data(mfile);
    molecule.add_implicit_hydrogens();
    molecule.generate_new_hydration();
    auto hist = molecule.get_histogram();

    auto make_fit = [&] () {
        std::shared_ptr<fitter::HydrationFitter> fitter;
        if (settings::hist::fit_excluded_volume) {
            fitter = std::make_shared<fitter::ExcludedVolumeFitter>(data, std::make_unique<hist::CompositeDistanceHistogramFFExplicit>(*static_cast<hist::CompositeDistanceHistogramFFExplicit*>(hist.get())));
        } else {
            fitter = std::make_shared<fitter::HydrationFitter>(data, std::make_unique<hist::CompositeDistanceHistogramFFExplicit>(*static_cast<hist::CompositeDistanceHistogramFFExplicit*>(hist.get())));
        }
        auto fit = fitter->fit();
        return fit->fval/fit->dof;
    };

    std::vector<double> chi2_voronoi;
    std::vector<double> chi2_mf;
    std::vector<constants::displaced_volume::detail::DisplacedVolumeSet> tables_voronoi;
    std::vector<constants::displaced_volume::detail::DisplacedVolumeSet> tables_mf;
    for (unsigned int i = 0; i < repeats; i++) {
        std::cout << 100*double(i)/repeats << "%\r" << std::flush;
        tables_voronoi  .push_back(generate_voronoi_exv_table());
        chi2_voronoi    .push_back(make_fit());
        tables_mf       .push_back(generate_mf_exv_table());
        chi2_mf         .push_back(make_fit());
    }
    std::cout << "Min voronoi is " << *std::min_element(chi2_voronoi.begin(), chi2_voronoi.end()) << std::endl;
    std::cout << "Min mf is " << *std::min_element(chi2_mf.begin(), chi2_mf.end()) << std::endl;

    std::ofstream out(settings::general::output + mfile.stem() + ".txt");
    out << "### Voronoi ###\n" << as_string(chi2_voronoi, tables_voronoi) << "\n\n";
    out << "### Minimum Fluctuation ###\n" << as_string(chi2_mf, tables_mf) << "\n\n";
    out.close();

    Axis axis(0, 1000, 10000);
    hist::Histogram chi2_hist_voronoi(axis);
    hist::Histogram chi2_hist_mf(axis);
    chi2_hist_voronoi.bin(chi2_voronoi);
    chi2_hist_mf.bin(chi2_mf);
    chi2_hist_voronoi.normalize();
    chi2_hist_mf.normalize();

    double count = 0, end = axis.bins;
    while (count < 0.01 && 0 < end) {
        count += chi2_hist_voronoi[--end];
        count += chi2_hist_mf[end];
    }

    if (100 < end) {
        unsigned int merge = end/100;
        chi2_hist_voronoi.merge(merge);
        chi2_hist_mf.merge(merge);
    }

    settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Voronoi_implicit_H;
    auto chi2_mean_voronoi = make_fit();
    std::stringstream ss; ss << "$\\chi^2_r = " << std::setprecision(1) << std::fixed << chi2_mean_voronoi << "$, Voronoi";
    std::string legend_voronoi = ss.str();

    settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::MinimumFluctutation_implicit_H;
    auto chi2_mean_mf = make_fit();
    ss.str(""); ss << "$\\chi^2_r = " << std::setprecision(1) << std::fixed << chi2_mean_mf << "$, Minimum Fluctuation";
    std::string legend_mf = ss.str();

    settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Traube;
    auto chi2_traube = make_fit();
    ss.str(""); ss << "$\\chi^2_r = " << std::setprecision(1) << std::fixed << chi2_traube << "$, Traube";
    std::string legend_traube = ss.str();

    end = std::max<int>(end, axis.get_bin(chi2_traube)+10);
    end = std::max<int>(end, axis.get_bin(chi2_mean_voronoi)+10);
    end = std::max<int>(end, axis.get_bin(chi2_mean_mf)+10);

    plots::PlotDataset()
        .plot(chi2_hist_voronoi.as_dataset(), plots::PlotOptions({{"legend", legend_voronoi}, {"marker", true}, {"lw", 0.5}, {"xlabel", "$\\chi^2_r$"}, {"ylabel", "density"}, {"xlim", Limit{0, axis.get_bin_value(end)}}, {"color", "tab:blue"}}))
        .plot(chi2_hist_mf.as_dataset(),      plots::PlotOptions({{"legend", legend_mf}, {"marker", true}, {"lw", 0.5}, {"color", "tab:orange"}}))
        .vline(chi2_mean_voronoi, plots::PlotOptions({{"color", "tab:blue"}, {"ls", style::line::dashed}}))
        .vline(chi2_mean_mf, plots::PlotOptions({{"color", "tab:orange"}, {"ls", style::line::dashed}}))
        .vline(chi2_traube, plots::PlotOptions({{"color", "r"}, {"legend", legend_traube}, {"ls", style::line::dashed}}))
        .save(settings::general::output + mfile.stem() + ".png");
    return 0;
}