#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distribution/WeightedDistribution.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <settings/MoleculeSettings.h>
#include <settings/GeneralSettings.h>

using namespace hist;
using namespace data;
using namespace data::record;

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

void set_unity_charge(data::Molecule& protein) {
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

std::vector<double> d_exact = {
    0, 
    std::sqrt(2), 
    2, 
    std::sqrt(8), 
    std::sqrt(12)
};

TEST_CASE("WeightedDistribution: tracks content") {
    hist::WeightedDistribution::reset();
    SECTION("simple") {
        hist::WeightedDistribution1D p(10, 0);
        auto width = constants::axes::d_axis.width();
        p.add(0, 1);
        p.add(width/2, 2);

        auto weighted_bins = WeightedDistribution::get_weighted_bins();
        REQUIRE(weighted_bins[0] == 0);
        REQUIRE(weighted_bins[1] == width/2);
        REQUIRE(weighted_bins[2] == 2*width);
    }

    SECTION("WeightedDistribution1D") {
        hist::WeightedDistribution1D p(10, 0);
        auto width = constants::axes::d_axis.width();
        p.add(0, 1);
        p.add(width/2, 2);
        p.add(3*width/4, 4);

        auto weighted_bins = WeightedDistribution::get_weighted_bins();
        REQUIRE(weighted_bins[0] == 0);
        REQUIRE(weighted_bins[1] == 2.5*width/4);
        REQUIRE(weighted_bins[2] == 2*width);
    }

    SECTION("WeightedDistribution2D") {
        hist::WeightedDistribution2D p(10, 10, 0);
        auto width = constants::axes::d_axis.width();
        p.add(0, 0, 1);
        p.add(1, width/2, 2);
        p.add(2, 3*width/4, 4);

        auto weighted_bins = WeightedDistribution::get_weighted_bins();
        REQUIRE(weighted_bins[0] == 0);
        REQUIRE(weighted_bins[1] == 2.5*width/4);
        REQUIRE(weighted_bins[2] == 2*width);
    }

    SECTION("WeightedDistribution3D") {
        hist::WeightedDistribution3D p(10, 10, 10, 0);
        auto width = constants::axes::d_axis.width();
        p.add(0, 0, 0, 1);
        p.add(1, 1, width/2, 2);
        p.add(2, 2, 3*width/4, 4);

        auto weighted_bins = WeightedDistribution::get_weighted_bins();
        REQUIRE(weighted_bins[0] == 0);
        REQUIRE(weighted_bins[1] == 2.5*width/4);
        REQUIRE(weighted_bins[2] == 2*width);
    }
}

TEST_CASE("CompositeDistanceHistogram::debye_transform (weighted)") {
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
                    16*std::sin(q_axis[q]*d_exact[1])/(q_axis[q]*d_exact[1]) +
                    24*std::sin(q_axis[q]*d_exact[2])/(q_axis[q]*d_exact[2]) + 
                    24*std::sin(q_axis[q]*d_exact[3])/(q_axis[q]*d_exact[3]) + 
                    8 *std::sin(q_axis[q]*d_exact[4])/(q_axis[q]*d_exact[4]);
                Iq_exp[q] += dsum*std::pow(ff(q_axis[q]), 2);
            }
        }

        {
            WeightedDistribution::reset();
            auto Iq = hist::HistogramManager<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            WeightedDistribution::reset();
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
            WeightedDistribution::reset();
            auto Iq = hist::HistogramManager<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
        {
            WeightedDistribution::reset();
            auto Iq = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
        }
    }
}
