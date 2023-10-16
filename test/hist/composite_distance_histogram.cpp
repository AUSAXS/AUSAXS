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
    Axis axis(1, 10, 1);
    return hist::CompositeDistanceHistogram(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(p), axis);
}

TEST_CASE("CompositeDistanceHistogram::reset_water_scaling_factor") {
    auto hist = generate_random(100);
    auto p = hist.get_total_counts();
    hist.apply_water_scaling_factor(2);
    CHECK(hist.get_total_counts() != p);
    hist.reset_water_scaling_factor();
    CHECK(hist.get_total_counts() == p);
}

TEST_CASE("CompositeDistanceHistogram::apply_water_scaling_factor") {
    settings::molecule::use_effective_charge = false;

    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    std::vector<Atom> b1 =   {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b2 =   {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b3 =   {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Water> w =   {Water(Vector3<double>(1, -1,  1), 1, constants::atom_t::C, "C", 1),  Water(Vector3<double>(1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3)};
    Molecule protein(a, w);

    auto hist = protein.get_histogram();
    std::vector<double> p_pp = hist->get_pp_counts();
    std::vector<double> p_hp = hist->get_hp_counts();
    std::vector<double> p_hh = hist->get_hh_counts();

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

TEST_CASE("CompositeDistanceHistogramFF::debye_transform") {
    settings::molecule::use_effective_charge = false;

    SECTION("no water") {
        std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b5 = {Atom(Vector3<double>( 0,  0,  0), 1, constants::atom_t::C, "C", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
        Molecule protein(a);

        set_unity_charge(protein);

        // check distances
        // calculation: 8 identical points. 
        //      each point has:
        //          1 line  of length 0
        //          3 lines of length 2
        //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
        //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
        // center point has:
        //          1 line  of length 0
        //          16 lines of length sqrt(2) (two for each point)

        std::vector<double> p_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};
        // std::vector<double> d =     {0, 2, std::sqrt(8), std::sqrt(12)};
        std::vector<double> d =     {0, 2.5, 2.5, 3.5, 1.5};
        std::vector<double> Iq_exp;
        {
            std::vector<double> q_axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
            Iq_exp.resize(q_axis.size(), 0);
            auto ff = [] (double q) {return std::exp(-q*q/2);};

            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                double dsum = 
                    9 + 
                    24*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
                    24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                    8*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) +
                    16*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

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

        // check distances
        // calculation: 8 identical points + 1 water in the middle. 
        //      each point has:
        //          1 line  of length 0
        //          1 line of length sqrt(2) = 1.41
        //          3 lines of length 2
        //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
        //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
        std::vector<double> p_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};
        // std::vector<double> d =     {0, 2, std::sqrt(8), std::sqrt(12), std::sqrt(2)};
        std::vector<double> d =     {0, 2.5, 2.5, 3.5, 1.5};
        std::vector<double> Iq_exp;
        {
            std::vector<double> q_axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
            Iq_exp.resize(q_axis.size(), 0);
            auto ff = [] (double q) {return std::exp(-q*q/2);};

            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                double aasum = 
                    8 + 
                    24*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
                    24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                    8*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]);
                Iq_exp[q] += aasum*std::pow(ff(q_axis[q]), 2);

                double awsum = 8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
                Iq_exp[q] += 2*awsum*std::pow(ff(q_axis[q]), 2);
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