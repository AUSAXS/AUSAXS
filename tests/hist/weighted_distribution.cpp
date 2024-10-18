#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <dataset/SimpleDataset.h>
#include <settings/MoleculeSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <grid/Grid.h>
#include <io/ExistingFile.h>
#include <table/ArrayDebyeTable.h>

#include "hist/hist_test_helper.h"
#include "hist/intensity_calculator/DistanceHistogram.h"
#include "plots/PlotDataset.h"

using namespace hist;
using namespace data;
using namespace data::record;

TEST_CASE("WeightedDistribution: tracks content") {
    SECTION("simple") {
        hist::WeightedDistribution1D p(10);
        auto width = constants::axes::d_axis.width();
        p.add(0, 1);
        p.add(width/2, 2);

        auto weighted_bins = p.get_weighted_axis();
        REQUIRE_THAT(weighted_bins[0], Catch::Matchers::WithinAbs(0, 1e-3));
        REQUIRE_THAT(weighted_bins[1], Catch::Matchers::WithinAbs(width/2, 1e-3));
        REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2*width, 1e-3));
    }

    SECTION("WeightedDistribution1D") {
        hist::WeightedDistribution1D p(10);
        auto width = constants::axes::d_axis.width();
        p.add(0, 1);
        p.add(width/2, 2);
        p.add(3*width/4, 4);

        auto weighted_bins = p.get_weighted_axis();
        REQUIRE_THAT(weighted_bins[0], Catch::Matchers::WithinAbs(0, 1e-3));
        REQUIRE_THAT(weighted_bins[1], Catch::Matchers::WithinAbs(2.5*width/4, 1e-3));
        REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2*width, 1e-3));
    }

    SECTION("WeightedDistribution2D") {
        hist::WeightedDistribution2D p(10, 10);
        auto width = constants::axes::d_axis.width();
        p.add(0, 0, 1);
        p.add(1, width/2, 2);
        p.add(2, 3*width/4, 4);

        auto weighted_bins = p.get_weighted_axis();
        REQUIRE_THAT(weighted_bins[0], Catch::Matchers::WithinAbs(0, 1e-3));
        REQUIRE_THAT(weighted_bins[1], Catch::Matchers::WithinAbs(2.5*width/4, 1e-3));
        REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2*width, 1e-3));
    }

    SECTION("WeightedDistribution3D") {
        hist::WeightedDistribution3D p(10, 10, 10);
        auto width = constants::axes::d_axis.width();
        p.add(0, 0, 0, 1);
        p.add(1, 1, width/2, 2);
        p.add(2, 2, 3*width/4, 4);

        auto weighted_bins = p.get_weights();
        REQUIRE_THAT(weighted_bins[0], Catch::Matchers::WithinAbs(0, 1e-3));
        REQUIRE_THAT(weighted_bins[1], Catch::Matchers::WithinAbs(2.5*width/4, 1e-3));
        REQUIRE_THAT(weighted_bins[2], Catch::Matchers::WithinAbs(2*width, 1e-3));
    }
}

class DistanceHistogramDebug : public DistanceHistogram {
    public:
        DistanceHistogramDebug(DistanceHistogram&& other) : DistanceHistogram(std::move(other)) {}
        auto get_sinc_table() const {return DistanceHistogram::get_sinc_table();}
};

TEST_CASE("WeightedDistribution: sinc_table") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b5 = {Atom(Vector3<double>( 0,  0,  0), 1, constants::atom_t::C, "C", 1)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
    Molecule protein(a);

    auto hist = hist::HistogramManagerMT<true>(&protein).calculate_all();
    auto Iq = hist->debye_transform();

    const auto& bins = constants::axes::d_vals;
    auto table = DistanceHistogramDebug(std::move(hist)).get_sinc_table();
    for (unsigned int q = 0; q < table->size_q(); ++q) {
        std::vector<double> sinc(20);
        for (unsigned int d = 0; d < 20; ++d) {
            double qd = constants::axes::q_vals[q]*bins[d];
            double val = 0;
            if (qd < 1e-3) {val = 1 - qd*qd/6 + qd*qd*qd*qd/120;}
            else {val = std::sin(qd)/qd;}
            REQUIRE_THAT(table->lookup(q, d), Catch::Matchers::WithinAbs(val, 1e-6));
            sinc[d] = val;
        }
        std::transform(sinc.begin(), sinc.end(), table->begin(q), sinc.begin(), std::minus<double>());
        REQUIRE_THAT(std::accumulate(sinc.begin(), sinc.end(), 0.0), Catch::Matchers::WithinAbs(0, 1e-6));
    }
}

// Check that the weighted distance axis is correctly calculated for all histogram managers.
TEST_CASE("WeightedDistribution: distance_calculators") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b5 = {Atom(Vector3<double>( 0,  0,  0), 1, constants::atom_t::C, "C", 1)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
    Molecule protein(a);

    { // hm
        CHECK(SimpleCube::check_default(hist::HistogramManager<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManager<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt
        CHECK(SimpleCube::check_default(hist::HistogramManagerMT<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMT<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt_ff_avg
        CHECK(SimpleCube::check_default(hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMTFFAvg<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt_ff_explicit
        CHECK(SimpleCube::check_default(hist::HistogramManagerMTFFExplicit<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMTFFExplicit<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt_ff_grid
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMTFFGrid(&protein).calculate_all()->get_d_axis()));
        // CHECK(check_exact(hist::HistogramManagerMTFFGrid<true>(&protein).calculate_all())); // exv cells dominates bin locs in this case
    }
    { // phm
        CHECK(SimpleCube::check_default(hist::PartialHistogramManager<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::PartialHistogramManager<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // phm_mt
        CHECK(SimpleCube::check_default(hist::PartialHistogramManagerMT<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::PartialHistogramManagerMT<true>(&protein).calculate_all()->get_d_axis()));
    }
}

// Check that the basic histogram managers agree on a weighted debye transform.
TEST_CASE("CompositeDistanceHistogram::debye_transform (weighted)") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::general::warnings = true;
    auto d_exact = SimpleCube::d_exact;

    SECTION("no water") {
        std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b5 = {Atom(Vector3<double>( 0,  0,  0), 1, constants::atom_t::C, "C", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
        Molecule protein(a);

        set_unity_charge(protein);

        std::vector<double> Iq_exp;
        {
            const auto& q_axis = constants::axes::q_vals;
            Iq_exp.resize(q_axis.size(), 0);
            auto ff2 = [] (double q) {return std::exp(-q*q);};

            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                double dsum = 
                    9 + 
                    16*std::sin(q_axis[q]*d_exact[1])/(q_axis[q]*d_exact[1]) +
                    24*std::sin(q_axis[q]*d_exact[2])/(q_axis[q]*d_exact[2]) + 
                    24*std::sin(q_axis[q]*d_exact[3])/(q_axis[q]*d_exact[3]) + 
                    8 *std::sin(q_axis[q]*d_exact[4])/(q_axis[q]*d_exact[4]);
                Iq_exp[q] += dsum*ff2(q_axis[q]);
            }
        }

        {
            auto Iq = hist::HistogramManager<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
    }

    SECTION("with water") {
        std::vector<Atom> b1 =  {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 =  {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 =  {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 =  {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Water> w = {Water(Vector3<double>( 0,  0,  0), 1, constants::atom_t::O, "HOH", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        DebugMolecule protein(a, w);

        set_unity_charge(protein);
        double Z = protein.get_volume_grid()*constants::charge::density::water/8;
        protein.set_volume_scaling(1./Z);

        std::vector<double> Iq_exp;
        {
            const auto& q_axis = constants::axes::q_vals;
            Iq_exp.resize(q_axis.size(), 0);
            auto ff = [] (double q) {return std::exp(-q*q/2);};

            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                double aasum = 
                    8 + 
                    24*std::sin(q_axis[q]*d_exact[2])/(q_axis[q]*d_exact[2]) + 
                    24*std::sin(q_axis[q]*d_exact[3])/(q_axis[q]*d_exact[3]) + 
                    8* std::sin(q_axis[q]*d_exact[4])/(q_axis[q]*d_exact[4]);
                Iq_exp[q] += aasum*std::pow(ff(q_axis[q]), 2);

                double awsum = 16*std::sin(q_axis[q]*d_exact[1])/(q_axis[q]*d_exact[1]);
                Iq_exp[q] += awsum*std::pow(ff(q_axis[q]), 2);
                Iq_exp[q] += 1*std::pow(ff(q_axis[q]), 2);
            }
        }

        {
            auto Iq = hist::HistogramManager<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
    }
}

#include <form_factor/ExvFormFactor.h>
#include <form_factor/FormFactor.h>
TEST_CASE("6lyz_exv", "[manual]") {
    auto exact = [] (const data::Molecule& molecule, double exv_radius) {
        container::Container2D<double> distances(molecule.get_atoms().size(), molecule.get_atoms().size());
        auto atoms = molecule.get_atoms();
        for (unsigned int i = 0; i < atoms.size(); ++i) {
            for (unsigned int j = 0; j < atoms.size(); ++j) {
                distances(i, j) = atoms[i].distance(atoms[j]);
            }
        }

        auto qaxis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        auto q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);
        form_factor::FormFactor ff = form_factor::ExvFormFactor(std::pow(2*exv_radius, 3));
        hist::ScatteringProfile I(qaxis);
        for (unsigned int q = q0; q < q0+qaxis.bins; ++q) {
            double sum = 0;
            for (unsigned int i = 0; i < atoms.size(); ++i) {
                for (unsigned int j = 0; j < atoms.size(); ++j) {
                    double qd = constants::axes::q_vals[q]*distances(i, j);
                    if (qd < 1e-6) {
                        // sum += std::pow(ff.evaluate(constants::axes::q_vals[q]), 2);
                        sum += 1;
                    } else {
                        // sum += std::pow(ff.evaluate(constants::axes::q_vals[q]), 2)*std::sin(qd)/(qd);
                        sum += std::sin(qd)/(qd);
                    }
                }
            }
            I.index(q-q0) = sum;
        }
        return I;
    };

    settings::molecule::use_effective_charge = false;
    settings::molecule::center = false;
    settings::axes::qmin = 5e-2;
    settings::axes::qmax = 1;
    settings::grid::exv::width = 1.5;
    settings::grid::min_exv_radius = 2;
    constants::radius::set_dummy_radius(2);

    data::Molecule protein("tests/files/6lyz_exv.pdb");
    for (auto& b : protein.get_bodies()) {for (auto& a : b.get_atoms()){a.element = constants::atom_t::dummy;}}
    auto Iq =  static_cast<hist::CompositeDistanceHistogramFFGrid*>(hist::HistogramManagerMTFFGrid(&protein).calculate_all().get())->get_profile_xx().as_dataset();
    auto Iqexact = exact(protein, settings::grid::exv::width).as_dataset();

    Iq.normalize();
    Iqexact.normalize();

    plots::PlotDataset()
        .plot(Iq, plots::PlotOptions(style::draw::line, {{"color", style::color::orange}, {"legend", "Unweighted"}, {"lw", 2}, {"yrange", Limit(1e-4, 1.1)}}))
        .plot(Iqexact, plots::PlotOptions(style::draw::line, {{"color", style::color::green}, {"legend", "Exact"}, {"ls", style::line::dashed}, {"lw", 2}}))
    .save("temp/tests/hist/6lyz_exv.png");
}

TEST_CASE("sphere_comparison", "[manual]") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::center = false;
    settings::axes::qmax = 1;
    auto lims = Limit3D(-50, 50, -50, 50, -50, 50);
    grid::Grid grid(lims);
    double radius = 15;
    double radius2 = radius*radius;
    auto axes = grid.get_axes();
    Vector3<double> center = grid.to_xyz(grid.get_center());
    for (unsigned int i = 0; i < axes.x.bins; ++i) {
        for (unsigned int j = 0; j < axes.y.bins; ++j) {
            for (unsigned int k = 0; k < axes.z.bins; ++k) {
                if (grid.to_xyz(i, j, k).distance2(center) < radius2) {
                    grid.grid.index(i, j, k) = grid::detail::VOLUME;
                }
            }
        }
    }
    auto loc = "temp/tests/hist/sphere.pdb";
    grid.save(loc);

    Molecule protein(loc);
    auto Iq1 = hist::HistogramManagerMT<false>(&protein).calculate_all()->debye_transform();
    auto Iq2 = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();

    plots::PlotDataset()
        .plot(Iq1.as_dataset(), plots::PlotOptions(style::draw::line, {{"color", style::color::orange}, {"legend", "Unweighted"}, {"lw", 2}}))
        .plot(Iq2.as_dataset(), plots::PlotOptions(style::draw::line, {{"color", style::color::blue}, {"legend", "Weighted"}, {"ls", style::line::dashed}, {"lw", 2}}))
    .save("temp/tests/hist/sphere_comparison.png");

    auto ratio = Iq1;
    for (unsigned int i = 0; i < ratio.get_counts().size(); ++i) {
        ratio.get_count(i) = Iq1.get_count(i)/Iq2.get_count(i);
    }

    plots::PlotDataset(
        ratio.as_dataset(),
        plots::PlotOptions(style::draw::line, {{"color", style::color::red}, {"legend", "Ratio"}, {"lw", 2}}))
    .save("temp/tests/hist/sphere_comparison_ratio.png");
}

TEST_CASE("real_comparison", "[manual]") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::center = false;
    settings::axes::qmax = 1;

    data::Molecule protein("tests/files/LAR1-2.pdb");
    auto Iq1 = hist::HistogramManagerMT<false>(&protein).calculate_all()->debye_transform();
    auto hist = hist::HistogramManagerMT<true>(&protein).calculate_all();
    auto Iq2 = hist->debye_transform();

    unsigned int counter = 0;
    auto weighted_bins = hist->get_d_axis();
    for (unsigned int i = 0; i < weighted_bins.size(); ++i) {
        if (weighted_bins[i] != constants::axes::d_vals[i]) {
            counter++;
        }
    }
    REQUIRE(counter != 0);

    {
        auto table = DistanceHistogramDebug(std::move(hist)).get_sinc_table();
        counter = 0;
        const auto& default_table = table::ArrayDebyeTable::get_default_table();
        REQUIRE(table->size_q() == default_table.size_q());
        REQUIRE(table->size_d() == default_table.size_d());
        for (unsigned int q = 0; q < table->size_q(); ++q) {
            for (unsigned int d = 0; d < table->size_d(); ++d) {
                if (table->lookup(q, d) != default_table.lookup(q, d)) {
                    counter++;
                }
            }
        }
        REQUIRE(counter != 0);
    }

    counter = 0;
    REQUIRE(Iq1.size() == Iq2.size());
    for (unsigned int i = 0; i < Iq1.size(); ++i) {
        if (Iq1.get_count(i) != Iq2.get_count(i)) {
            counter++;
        }
    }
    REQUIRE(counter != 0);

    plots::PlotDataset()
        .plot(Iq1.as_dataset(), plots::PlotOptions(style::draw::line, {{"color", style::color::orange}, {"legend", "Unweighted"}, {"lw", 2}}))
        .plot(Iq2.as_dataset(), plots::PlotOptions(style::draw::line, {{"color", style::color::blue}, {"legend", "Weighted"}, {"ls", style::line::dashed}, {"lw", 2}}))
    .save("temp/tests/hist/real_comparison.png");
}