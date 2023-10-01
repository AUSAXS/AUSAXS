#include <catch2/catch_test_macros.hpp>

#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/HistogramManagerMT.h>
#include <hist/HistogramManagerMTFF.h>
#include <hist/Histogram.h>
#include <hist/detail/FormFactor.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>
#include <hist/DebyeLookupTable.h>
#include <settings/All.h>

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

void set_unity_charge(Protein& protein) {
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

#define DEBYE_DEBUG 0
TEST_CASE("CompositeDistanceHistogramFF::debye_transform") {
    settings::protein::use_effective_charge = false;
    auto ff_carbon = hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::NEUTRAL_CARBON);
    auto ff_exv = hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::EXCLUDED_VOLUME);
    auto ff_w = hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::NEUTRAL_OXYGEN);
    std::vector<double> q_axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
    std::vector<double> Iq_exp(q_axis.size(), 0);

    SECTION("no water") {
        std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b5 = {Atom(Vector3<double>( 0,  0,  0), 1, constants::atom_t::C, "C", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
        Protein protein(a);

        set_unity_charge(protein);
        double Z = protein.get_excluded_volume()*constants::charge::density::water/9;
        protein.set_excluded_volume_scaling(1./Z);

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
        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double dsum = 
                9 + 
                24*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                8*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) +
                16*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
            Iq_exp[q] += dsum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);
            Iq_exp[q] -= 2*dsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);
            Iq_exp[q] += dsum*std::pow(ff_exv.evaluate(q_axis[q]), 2);
        }

        auto Iq = hist::HistogramManagerMTFF(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("with water") {
        std::vector<Atom> b1 =  {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 =  {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 =  {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 =  {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Water> w = {Water(Vector3<double>( 0,  0,  0), 1, constants::atom_t::O, "HOH", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Protein protein(a, w);

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
        unsigned int qcheck = 1;
        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double aasum = 
                8 + 
                24*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                8*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]);

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "8*1" << std::endl;
                    std::cout << "8*sin(" << q_axis[q] << "*" << d[4] << ")/(" << q_axis[q] << "*" << d[4] << ") = " << std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]) << std::endl;
                    std::cout << "24*sin(" << q_axis[q] << "*" << d[1] << ")/(" << q_axis[q] << "*" << d[1] << ") = " << std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) << std::endl;
                    std::cout << "24*sin(" << q_axis[q] << "*" << d[2] << ")/(" << q_axis[q] << "*" << d[2] << ") = " << std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) << std::endl;
                    std::cout << "8*sin(" << q_axis[q] << "*" << d[3] << ")/(" << q_axis[q] << "*" << d[3] << ") = " << std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) << std::endl;
                }
            #endif
            Iq_exp[q] += aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);                  // + aa
            Iq_exp[q] -= 2*aasum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);  // -2ax
            Iq_exp[q] += aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2);                     // + xx

            double awsum = 8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
            Iq_exp[q] += 2*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);    // +2aw
            Iq_exp[q] -= 2*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);       // -2ew
            Iq_exp[q] += 1*std::pow(ff_w.evaluate(q_axis[q]), 2);                           // + ww

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "aasum = " << aasum << std::endl;
                    std::cout << "awsum = " << awsum << std::endl;
                    std::cout << "(aa) Iq_exp[" << q << "] += " << aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(ae) Iq_exp[" << q << "] -= " << 2*aasum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(aw) Iq_exp[" << q << "] += " << 2*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(ee) Iq_exp[" << q << "] += " << aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(ew) Iq_exp[" << q << "] -= " << 2*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(ww) Iq_exp[" << q << "] += " << 1*std::pow(ff_w.evaluate(q_axis[q]), 2) << std::endl;
                }
            #endif
        }
        auto Iq = hist::HistogramManagerMTFF(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("real scalings") {
        std::vector<Atom> b1 =  {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 =  {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 =  {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 =  {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Water> w = {Water(Vector3<double>( 0,  0,  0), 1, constants::atom_t::O, "HOH", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Protein protein(a, w);

        double ZX = protein.get_excluded_volume()*constants::charge::density::water/8;
        double ZC = 6; // 6
        double ZO = 8; // 8

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
        unsigned int qcheck = 1;
        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double aasum = 
                8 + 
                24*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                8*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]);
            Iq_exp[q] += ZC*ZC*aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);                    // + aa
            Iq_exp[q] -= 2*ZC*ZX*aasum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);    // -2ax
            Iq_exp[q] += ZX*ZX*aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2);                       // + xx

            double awsum = 8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
            Iq_exp[q] += 2*ZC*ZO*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);      // +2aw
            Iq_exp[q] -= 2*ZO*ZX*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);         // -2ew
            Iq_exp[q] += 1*ZO*ZO*std::pow(ff_w.evaluate(q_axis[q]), 2);                             // + ww
        }
        auto Iq = hist::HistogramManagerMTFF(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("real data") {
        Protein protein("test/files/2epe.pdb");
        double ZX = protein.get_excluded_volume()*constants::charge::density::water/protein.atom_size();
        double ZC = 6;
        double ZO = 8;
        for (auto& body : protein.get_bodies()) {
            for (auto& atom : body.get_atoms()) {
                atom.set_element(constants::atom_t::C);
                atom.set_effective_charge(ZC);
            }
        }
        protein.generate_new_hydration();
        for (auto& water : protein.get_waters()) {
            water.set_effective_charge(ZO);
        }
        auto Iq = hist::HistogramManagerMTFF(&protein).calculate_all()->debye_transform();

        unsigned int N = protein.atom_size();
        unsigned int M = protein.water_size();
        REQUIRE(utility::approx(Iq[0], std::pow(N*ZC + M*ZO - N*ZX, 2)));
    }
}