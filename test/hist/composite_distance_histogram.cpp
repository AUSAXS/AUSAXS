#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>
#include <constants/Constants.h>
#include <settings/All.h>

#include "../test/hist/hist_test_helper.h"

using namespace data::record;
using namespace data;


hist::CompositeDistanceHistogram generate_random(unsigned int size) {
    hist::Distribution1D p_pp(size), p_hp(size), p_hh(size), p(size);
    for (unsigned int i = 0; i < size; ++i) {
        p_pp.index(i) = rand() % 100;
        p_hp.index(i) = rand() % 100;
        p_hh.index(i) = rand() % 100;
        p.index(i) = p_pp.index(i) + 2*p_hp.index(i) + p_hh.index(i);
    }
    Axis axis(1, 10, size);
    return hist::CompositeDistanceHistogram(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(p), axis);
}

TEST_CASE("CompositeDistanceHistogram::reset_water_scaling_factor") {
    settings::general::warnings = false;
    auto hist = generate_random(100);
    auto p = hist.get_total_counts();
    hist.apply_water_scaling_factor(2);
    CHECK(hist.get_total_counts() != p);
    hist.reset_water_scaling_factor();
    CHECK(hist.get_total_counts() == p);
}

TEST_CASE("CompositeDistanceHistogram::apply_water_scaling_factor") {
    settings::general::warnings = false;
    settings::molecule::use_effective_charge = false;

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

    hist->reset_water_scaling_factor();
    for (unsigned int i = 0; i < p_pp.size(); i++) {
        REQUIRE_THAT(p_pp[i] + 2*p_hp[i] + p_hh[i], Catch::Matchers::WithinRel(hist->get_total_counts()[i]));
    }
}

TEST_CASE("CompositeDistanceHistogram::debye_transform") {
    settings::molecule::use_effective_charge = false;
    settings::general::warnings = true;

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
        }

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
                    24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                    24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                    8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
                Iq_exp[q] += aasum*std::pow(ff(q_axis[q]), 2);

                double awsum = 16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]);
                Iq_exp[q] += awsum*std::pow(ff(q_axis[q]), 2);
                Iq_exp[q] += 1*std::pow(ff(q_axis[q]), 2);
            }
        }

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
}

TEST_CASE("qmin & qmax") {
    Molecule protein("test/files/2epe.pdb");
    {
        auto Iqf = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform();
        settings::axes::qmin = 1e-2;
        settings::axes::qmax = 0.7;
        auto Iqr = hist::HistogramManager<false>(&protein).calculate_all()->debye_transform();

        auto start = Iqf.get_axis().get_bin(settings::axes::qmin);
        auto end = Iqf.get_axis().get_bin(settings::axes::qmax);
        std::vector<double> Iqf_r(Iqf.get_counts().begin()+start, Iqf.get_counts().begin()+end+1);
        REQUIRE(compare_hist(Iqf_r, Iqr));
    }
}

TEST_CASE("CompositeDistanceHistogram::get_profile") {
    data::Molecule protein("test/files/2epe.pdb");
    auto hist = hist::HistogramManager<false>(&protein).calculate_all();
    auto Iq = hist->debye_transform();
    auto profile_sum = hist->get_profile_aa() + hist->get_profile_aw() + hist->get_profile_ww();
    REQUIRE(Iq.size() == profile_sum.size());
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        REQUIRE_THAT(Iq[i], Catch::Matchers::WithinAbs(profile_sum[i], 1e-6));
    }
}