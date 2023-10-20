#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/CompositeDistanceHistogram.h>
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

using namespace data::record;
using namespace data;

hist::CompositeDistanceHistogram generate_random(unsigned int size) {
    std::vector<double> p_pp(size), p_hp(size), p_hh(size), p(size);
    for (unsigned int i = 0; i < size; ++i) {
        p_pp[i] = rand() % 100;
        p_hp[i] = rand() % 100;
        p_hh[i] = rand() % 100;
        p[i] = p_pp[i] + 2*p_hp[i] + p_hh[i];
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

void set_unity_charge(Molecule& protein) {
    // set the weights to 1 so we can analytically determine the result
    // waters
    for (auto& atom : protein.get_waters()) {
        atom.set_effective_charge(1);
    }
    // atoms
    for (auto& body : protein.get_bodies()) {
        for (auto& atom : body.get_atoms()) {
            atom.set_effective_charge(1);
        }
    }
}

bool compare_hist(Vector<double> p1, Vector<double> p2) {
    unsigned int pmin = std::min(p1.size(), p2.size());
    for (unsigned int i = 0; i < pmin; i++) {
        if (!utility::approx(p1[i], p2[i], 1e-6, 0.01)) {
            std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
            return false;
        }
    }

    return true;
}

// calculation: 8 points
//          1 line  of length 0
//          3 lines of length 2
//          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
//          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
//
// calculation: 1 center point
//          1 line  of length 0
//          16 lines of length sqrt(2) = 1.41 (counting both directions)
//
// sum:
//          9 line  of length 0
//          16 lines of length sqrt(2)
//          24 lines of length 2
//          24 lines of length sqrt(8)
//          8 lines of length sqrt(12)
auto width = constants::axes::d_axis.width();
std::vector<double> d = {
    0, 
    constants::axes::d_vals[std::round(std::sqrt(2)/width)], 
    constants::axes::d_vals[std::round(2./width)], 
    constants::axes::d_vals[std::round(std::sqrt(8)/width)], 
    constants::axes::d_vals[std::round(std::sqrt(12)/width)]
};

TEST_CASE("CompositeDistanceHistogramFF::debye_transform") {
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
            auto Iq = hist::HistogramManager(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::HistogramManagerMT(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::PartialHistogramManager(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::PartialHistogramManagerMT(&protein).calculate_all()->debye_transform();
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
        Molecule protein(a, w);

        set_unity_charge(protein);
        double Z = protein.get_excluded_volume()*constants::charge::density::water/8;
        protein.set_excluded_volume_scaling(1./Z);

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
            auto Iq = hist::HistogramManager(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::HistogramManagerMT(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::PartialHistogramManager(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            auto Iq = hist::PartialHistogramManagerMT(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
    }
}

TEST_CASE("CompositeDistanceHistogram::get_profile") {
    data::Molecule protein("test/files/2epe.pdb");
    auto hist = hist::HistogramManager(&protein).calculate_all();
    auto Iq = hist->debye_transform();
    auto profile_sum = hist->get_profile_aa() + hist->get_profile_aw() + hist->get_profile_ww();
    REQUIRE(Iq.size() == profile_sum.size());
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        REQUIRE_THAT(Iq[i], Catch::Matchers::WithinAbs(profile_sum[i], 1e-6));
    }
}