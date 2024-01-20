#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <hist/distribution/WeightedDistribution.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/HistFwd.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <settings/MoleculeSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <hydrate/Grid.h>
#include <io/ExistingFile.h>
#include <plots/All.h>
#include <table/ArrayDebyeTable.h>
#include <table/VectorDebyeTable.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/FormFactor.h>
#include <utility/Utility.h>

#include <sstream>
#include <filesystem>
#include <ranges>
#include <cassert>

auto get_profile_xx = [] (hist::CompositeDistanceHistogramFFGrid* h) {
    auto weighted_bins = hist::WeightedDistribution::get_weighted_bins();
    auto sinqd_table = table::VectorDebyeTable(weighted_bins);
    auto cp_aa = h->get_aa_counts_ff();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double xx_sum = std::inner_product(cp_aa.begin(form_factor::exv_bin, form_factor::exv_bin), cp_aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table.begin(q), 0.0);
        Iq[q-q0] += xx_sum;
    }
    return hist::ScatteringProfile(Iq, debye_axis);
};

auto get_profile_aa = [] (hist::CompositeDistanceHistogramFFGrid* h) {
    auto weighted_bins = hist::WeightedDistribution::get_weighted_bins();
    auto sinqd_table = table::VectorDebyeTable(weighted_bins);
    auto cp_aa = h->get_aa_counts_ff();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < form_factor::get_count_without_excluded_volume(); ++j) {
                double aa_sum = std::inner_product(cp_aa.begin(i, j), cp_aa.end(i, j), sinqd_table.begin(q), 0.0);
                Iq[q-q0] += aa_sum;
            }
        }
    }
    return hist::ScatteringProfile(Iq, debye_axis);
};

auto exact = [] (const data::Molecule& molecule) {
    container::Container2D<double> distances(molecule.get_atoms().size(), molecule.get_atoms().size());
    auto atoms = molecule.get_atoms();
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        for (unsigned int j = 0; j < atoms.size(); ++j) {
            distances(i, j) = atoms[i].distance(atoms[j]);
        }
    }

    auto qaxis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    auto q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
    hist::ScatteringProfile I(qaxis);
    for (unsigned int q = q0; q < q0+qaxis.bins; ++q) {
        double sum = 0;
        for (unsigned int i = 0; i < atoms.size(); ++i) {
            for (unsigned int j = 0; j < atoms.size(); ++j) {
                double qd = constants::axes::q_vals[q]*distances(i, j);
                if (qd < 1e-6) {
                    sum += atoms[i].effective_charge*atoms[j].effective_charge;
                } else {
                    sum += atoms[i].effective_charge*atoms[j].effective_charge*std::sin(qd)/(qd);
                }
            }
        }
        I.index(q-q0) = sum;
    }
    return I;
};

// Instructions: go vary constants::axes::d_axis bins a couple of times while running this program inbetween. 
// When enough data has been collected, call it with the argument "plot" to plot the results.
int main(int argc, char const *argv[]) {
    bool setup = false;
    bool structured = true;
    if (argc > 1 && std::string(argv[1]) == "setup") {
        setup = true;
    } else if (argc > 1 && std::string(argv[1]) == "plot") {
        structured = false;
    } else if (argc > 1 && std::string(argv[1]) == "plotx") {
        structured = true;
    } else {
        std::cout << "Usage: " << argv[0] << " [setup|plot|plotx]" << std::endl;
        return 1;
    }

    settings::molecule::use_effective_charge = false;
    settings::molecule::center = false;
    settings::axes::qmin = 5e-2;
    settings::axes::qmax = 1;
    settings::grid::exv_radius = 1.5;
    settings::grid::rvol = 2;
    constants::radius::set_dummy_radius(2);
    if (setup) {
        data::Molecule protein_6lyz_exv("test/files/6lyz_exv.pdb");
        for (auto& b : protein_6lyz_exv.get_bodies()) {for (auto& a : b.get_atoms()){a.element = constants::atom_t::dummy;}}
        auto Iq =  get_profile_xx(static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid<false>(&protein_6lyz_exv).calculate_all().get())).as_dataset();
        auto Iqw = get_profile_xx(static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid<true>(&protein_6lyz_exv).calculate_all().get())).as_dataset();

        hist::WeightedDistribution::reset();
        data::Molecule protein_6lyz("test/files/6lyz.pdb");
        auto Iq2 =  get_profile_aa(static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid<false>(&protein_6lyz).calculate_all().get())).as_dataset();
        auto Iqw2 = get_profile_aa(static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid<true>(&protein_6lyz).calculate_all().get())).as_dataset();

        Iq.normalize();
        Iqw.normalize();
        Iq2.normalize();
        Iqw2.normalize();

        std::string width;
        {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << constants::axes::d_axis.width();
            width = ss.str();
        }

        Iq.save(  "output/bin_size_analysis/files/unweighted_structured_"   + width + ".dat");
        Iqw.save( "output/bin_size_analysis/files/weighted_structured_"     + width + ".dat");
        Iq2.save( "output/bin_size_analysis/files/unweighted_unstructured_" + width + ".dat");
        Iqw2.save("output/bin_size_analysis/files/weighted_unstructured_"   + width + ".dat");
    } else {
        data::Molecule protein(structured ? "test/files/6lyz_exv.pdb" : "test/files/6lyz.pdb");
        if (structured) {for (auto& b : protein.get_bodies()) {for (auto& a : b.get_atoms()){a.set_effective_charge(1);}}}
        protein.clear_hydration();
        auto Iqexact = exact(protein).as_dataset();
        Iqexact.normalize();

        settings::axes::qmin = 1e-4;

        std::vector<std::string> files_unweighted;
        std::vector<std::string> files_weighted;
        std::string prefix = structured ? "_structured" : "_unstructured";
        for (auto& p : std::filesystem::directory_iterator("output/bin_size_analysis/files")) {
            if (p.path().extension() == ".dat") {
                if (p.path().stem().string().find("unweighted" + prefix) != std::string::npos) {
                    files_unweighted.push_back(p.path().string());
                } else if (p.path().stem().string().find("weighted" + prefix) != std::string::npos) {
                    files_weighted.push_back(p.path().string());
                }
            }
        }
        std::sort(files_unweighted.begin(), files_unweighted.end(), std::greater<>());
        std::sort(files_weighted.begin(), files_weighted.end(), std::greater<>());

        std::vector<SimpleDataset> weighted;
        for (auto& f : files_weighted) {
            weighted.push_back(SimpleDataset(f));
        }

        std::vector<SimpleDataset> unweighted;
        for (auto& f : files_unweighted) {
            unweighted.push_back(SimpleDataset(f));
        }

        plots::PlotDataset p_all;
        plots::PlotDataset p_weighted, p_unweighted;
        plots::PlotDataset p_intensity_weighted;
        p_intensity_weighted.plot(Iqexact, plots::PlotOptions({{"color", style::color::black}, {"legend", "exact"}, {"lw", 2}, {"ls", style::line::solid}, {"logy", true}, {"logx", true}}));
        for (unsigned int i = 0; i < weighted.size(); ++i) {
            auto color = style::color::next();
            std::string width;
            std::ranges::copy(files_unweighted[i] | std::views::filter([] (char c) {return std::isdigit(c) || c == '.';}), std::back_inserter(width));
            width.resize(width.size()-1); // remove the last '.' character

            p_intensity_weighted.plot(weighted[i], plots::PlotOptions({{"color", color}, {"legend", "width " + width}, {"lw", 2}, {"ls", style::line::dashed}, {"logy", true}, {"logx", true}}));
            for (unsigned int j = 0; j < weighted[i].size(); ++j) {
                assert(utility::approx(  weighted[i].x(j), Iqexact.x(j), 1e-6, 0));
                assert(utility::approx(unweighted[i].x(j), Iqexact.x(j), 1e-6, 0));
                weighted[i].y(j) /= std::abs(Iqexact.y(j));
                unweighted[i].y(j) /= std::abs(Iqexact.y(j));
            }

            p_unweighted.plot(unweighted[i], plots::PlotOptions({{"color", color}, {"legend", "width " + width}, {"lw", 2}, {"yrange", Limit(0, 2)}}));
            p_weighted.plot(    weighted[i], plots::PlotOptions({{"color", color}, {"legend", "width " + width}, {"lw", 2}, {"yrange", Limit(0, 2)}}));
            p_all.plot(       unweighted[i], plots::PlotOptions({{"color", color}, {"legend", "unweighted, width " + width}, {"lw", 2}, {"ls", style::line::dashed}}));
            p_all.plot(         weighted[i], plots::PlotOptions({{"color", color}, {"legend", "weighted, width " + width}, {"lw", 2}, {"yrange", Limit(0, 2)}}));
        }
        p_all.hline(1, plots::PlotOptions({{"color", style::color::black}, {"lw", 1}}));
        p_weighted.hline(1, plots::PlotOptions({{"color", style::color::black}, {"lw", 1}}));
        p_unweighted.hline(1, plots::PlotOptions({{"color", style::color::black}, {"lw", 1}}));

        p_all.save("output/bin_size_analysis/all" + prefix + ".png");
        p_weighted.save("output/bin_size_analysis/weighted" + prefix + ".png");
        p_unweighted.save("output/bin_size_analysis/unweighted" + prefix + ".png");
        p_intensity_weighted.save("output/bin_size_analysis/intensity_weighted" + prefix + ".png");
    }
}