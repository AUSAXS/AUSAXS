#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <form_factor/FormFactor.h>
#include <dataset/SimpleDataset.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>
#include <constants/Constants.h>
#include <settings/All.h>
#include <plots/All.h>

#include "hist/hist_test_helper.h"

using namespace data::record;
using namespace data;

// calculate the exact aa profile assuming all atoms are carbon
auto exact_aa_carbon = [] (const data::Molecule& molecule) {
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
                    sum += 1;
                } else {
                    sum += std::sin(qd)/(qd);
                }
            }
        }
        I.index(q-q0) = sum;
    }
    return I;
};

hist::CompositeDistanceHistogram generate_random(unsigned int size) {
    hist::Distribution1D p_pp(size), p_hp(size), p_hh(size), p(size);
    for (unsigned int i = 0; i < size; ++i) {
        p_pp.index(i) = rand() % 100;
        p_hp.index(i) = rand() % 100;
        p_hh.index(i) = rand() % 100;
        p.index(i) = p_pp.index(i) + 2*p_hp.index(i) + p_hh.index(i);
    }
    Axis axis(1, 10, size);
    return hist::CompositeDistanceHistogram(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(p));
}

TEST_CASE("CompositeDistanceHistogram::reset_water_scaling_factor") {
    settings::general::warnings = false;
    auto hist = generate_random(100);
    auto p = hist.get_total_counts();
    hist.apply_water_scaling_factor(2);
    CHECK(hist.get_total_counts() != p);
    hist.apply_water_scaling_factor(1);
    CHECK(hist.get_total_counts() == p);
}

TEST_CASE("CompositeDistanceHistogram::apply_water_scaling_factor") {
    settings::general::warnings = false;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;

    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    std::vector<Atom> b1 =   {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b2 =   {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b3 =   {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Water> w =   {Water(Vector3<double>(1, -1,  1), 1, constants::atom_t::C, "C", 1),  Water(Vector3<double>(1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3)};
    Molecule protein(a, w);

    auto hist = protein.get_histogram();
    std::vector<double> p_pp = hist->get_aa_counts();
    std::vector<double> p_hp = hist->get_aw_counts();
    std::vector<double> p_hh = hist->get_ww_counts();

    hist->apply_water_scaling_factor(2);
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 2*2*p_hp[i] + 4*p_hh[i], Catch::Matchers::WithinRel(hist->get_total_counts()[i]));
    }

    hist->apply_water_scaling_factor(3);
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 3*2*p_hp[i] + 9*p_hh[i], Catch::Matchers::WithinRel(hist->get_total_counts()[i]));
    }

    hist->apply_water_scaling_factor(1);
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 2*p_hp[i] + p_hh[i], Catch::Matchers::WithinRel(hist->get_total_counts()[i]));
    }
}

TEST_CASE("CompositeDistanceHistogram::debye_transform", "[files]") {
    settings::general::warnings = true;
    settings::general::verbose = false;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    auto d = SimpleCube::d;

    SECTION("no water") {
        std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b5 = {Atom(Vector3<double>( 0,  0,  0), 1, constants::atom_t::C, "C", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
        Molecule protein(a);

        set_unity_charge(protein);

        auto test_func = [&] (const auto& q_axis) {
            std::vector<double> Iq_exp;
            Iq_exp.resize(q_axis.size(), 0);
            auto ff = [] (double q) {return std::exp(-q*q/2);};

            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                double dsum = 
                    9 + 
                    16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) +
                    24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                    24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                    8 *std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
                Iq_exp[q] += dsum*std::pow(ff(q_axis[q]), 2);
            }
            return Iq_exp;
        };

        SECTION("default q-axis") {
            auto Iq_exp = test_func(constants::axes::q_vals);
            {
                auto Iq = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
            {
                auto Iq = hist::HistogramManagerMT<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
            {
                auto Iq = hist::PartialHistogramManager<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
            {
                auto Iq = hist::PartialHistogramManagerMT<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
        }

        SECTION("custom q-axis") {
            std::vector<double> q_axis(100);
            for (unsigned int i = 0; i < q_axis.size(); ++i) {
                q_axis[i] = (i+1)*0.1;
            }
            auto Iq_exp = test_func(q_axis);

            {
                auto Iq = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
            {
                auto Iq = hist::HistogramManagerMT<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
            {
                auto Iq = hist::PartialHistogramManager<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
            {
                auto Iq = hist::PartialHistogramManagerMT<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
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

        auto test_func = [&] (const auto& q_axis) {
            std::vector<double> Iq_exp;
            Iq_exp.resize(q_axis.size(), 0);
            auto ff = [] (double q) {return std::exp(-q*q/2);};

            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                double aasum = 
                    8 + 
                    24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                    24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                    8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
                Iq_exp[q] += aasum*std::pow(ff(q_axis[q]), 2);

                double awsum = 16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]);
                Iq_exp[q] += awsum*std::pow(ff(q_axis[q]), 2);
                Iq_exp[q] += 1*std::pow(ff(q_axis[q]), 2);
            }
            return Iq_exp;
        };

        SECTION("default q-axis") {
            auto Iq_exp = test_func(constants::axes::q_vals);
            {
                auto Iq = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
            {
                auto Iq = hist::HistogramManagerMT<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
            {
                auto Iq = hist::PartialHistogramManager<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
            {
                auto Iq = hist::PartialHistogramManagerMT<false>(&protein).calculate_all()->debye_transform();
                REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
            }
        }

        SECTION("custom q-axis") {
            std::vector<double> q_axis(100);
            for (unsigned int i = 0; i < q_axis.size(); ++i) {
                q_axis[i] = (i+1)*0.1;
            }
            auto Iq_exp = test_func(q_axis);

            {
                auto Iq = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
            {
                auto Iq = hist::HistogramManagerMT<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
            {
                auto Iq = hist::PartialHistogramManager<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
            {
                auto Iq = hist::PartialHistogramManagerMT<false>(&protein).calculate_all()->debye_transform(q_axis);
                REQUIRE(compare_hist(Iq_exp, Iq.y()));
            }
        }
    }

    SECTION("analytical") {
        auto data = GENERATE(
            "2epe",
            "6lyz",
            "c60",
            "diamond",
            "LAR1-2"
        );

        SECTION("real data (" + std::string(data) + ")") {
            data::Molecule protein("tests/files/" + std::string(data) + ".pdb");
            protein.clear_hydration();
            for (auto& body : protein.get_bodies()) {
                for (auto& atom : body.get_atoms()) {
                    atom.set_element(constants::atom_t::C);
                    atom.set_residue_name("UNK");
                    atom.set_group_name("C");
                    atom.set_effective_charge(1);
                }
            }

            auto exact = exact_aa_carbon(protein);
            auto ff = form_factor::storage::atomic::C;
            { // hm
                auto hm = hist::HistogramManager<true>(&protein).calculate_all()->get_profile_aa();
                auto axis = hm.get_axis().as_vector();
                std::transform(hm.get_counts().begin(), hm.get_counts().end(), axis.begin(), hm.get_counts().begin(), [] (double x, double q) {return x*std::exp(q*q);});
                REQUIRE(compare_hist(exact, hm, 0, 1e-2)); // 1% error allowed
            }
            { // hm_mt
                auto hm_mt = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_profile_aa();
                auto axis = hm_mt.get_axis().as_vector();
                std::transform(hm_mt.get_counts().begin(), hm_mt.get_counts().end(), axis.begin(), hm_mt.get_counts().begin(), [] (double x, double q) {return x*std::exp(q*q);});
                REQUIRE(compare_hist(exact, hm_mt, 0, 1e-2));
            }
            { // hm_mt_ff_avg
                auto hm_mt_ff_avg = hist::HistogramManagerMTFFAvg<true>(&protein).calculate_all()->get_profile_aa();
                auto axis = hm_mt_ff_avg.get_axis().as_vector();
                std::transform(hm_mt_ff_avg.get_counts().begin(), hm_mt_ff_avg.get_counts().end(), axis.begin(), hm_mt_ff_avg.get_counts().begin(), 
                    [ff] (double x, double q) {return x/std::pow(ff.evaluate(q), 2);}
                );
                REQUIRE(compare_hist(exact, hm_mt_ff_avg, 0, 1e-2));
            }
            { // hm_mt_ff_explicit
                auto hm_mt_ff_explicit = hist::HistogramManagerMTFFExplicit<true>(&protein).calculate_all()->get_profile_aa();
                auto axis = hm_mt_ff_explicit.get_axis().as_vector();
                std::transform(hm_mt_ff_explicit.get_counts().begin(), hm_mt_ff_explicit.get_counts().end(), axis.begin(), hm_mt_ff_explicit.get_counts().begin(), 
                    [ff] (double x, double q) {return x/std::pow(ff.evaluate(q), 2);}
                );
                REQUIRE(compare_hist(exact, hm_mt_ff_explicit, 0, 1e-2));
            }
            { // hm_mt_ff_grid
                auto hm_mt_ff_grid = hist::HistogramManagerMTFFGrid(&protein).calculate_all()->get_profile_aa();
                auto axis = hm_mt_ff_grid.get_axis().as_vector();
                std::transform(hm_mt_ff_grid.get_counts().begin(), hm_mt_ff_grid.get_counts().end(), axis.begin(), hm_mt_ff_grid.get_counts().begin(), 
                    [ff] (double x, double q) {return x/std::pow(ff.evaluate(q), 2);}
                );
                REQUIRE(compare_hist(exact, hm_mt_ff_grid, 0, 1e-2));
            }
            { // phm
                auto phm = hist::PartialHistogramManager<true>(&protein).calculate_all()->get_profile_aa();
                auto axis = phm.get_axis().as_vector();
                std::transform(phm.get_counts().begin(), phm.get_counts().end(), axis.begin(), phm.get_counts().begin(), [] (double x, double q) {return x*std::exp(q*q);});
                REQUIRE(compare_hist(exact, phm, 0, 1e-2));
            }
            { // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT<true>(&protein).calculate_all()->get_profile_aa();
                auto axis = phm_mt.get_axis().as_vector();
                std::transform(phm_mt.get_counts().begin(), phm_mt.get_counts().end(), axis.begin(), phm_mt.get_counts().begin(), [] (double x, double q) {return x*std::exp(q*q);});
                REQUIRE(compare_hist(exact, phm_mt, 0, 1e-2));
            }
        }
    }
}

TEST_CASE("qmin & qmax") {
    settings::general::verbose = false;
    Molecule protein("tests/files/2epe.pdb");
    auto Iqf = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform();
    settings::axes::qmin = 1e-2;
    settings::axes::qmax = 0.7;
    auto Iqr = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform();

    auto start = Iqf.get_axis().get_bin(settings::axes::qmin);
    auto end = Iqf.get_axis().get_bin(settings::axes::qmax);
    std::vector<double> Iqf_r(Iqf.get_counts().begin()+start, Iqf.get_counts().begin()+end);
    REQUIRE(compare_hist(Iqf_r, Iqr));
}

TEST_CASE("CompositeDistanceHistogram::get_profile") {
    settings::general::verbose = false;
    data::Molecule protein("tests/files/2epe.pdb");
    auto hist = hist::HistogramManager<false>(&protein).calculate_all();
    auto Iq = hist->debye_transform();
    auto profile_sum = hist->get_profile_aa() + hist->get_profile_aw() + hist->get_profile_ww();
    REQUIRE(Iq.size() == profile_sum.size());
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        REQUIRE_THAT(Iq[i], Catch::Matchers::WithinAbs(profile_sum[i], 1e-6));
    }
}