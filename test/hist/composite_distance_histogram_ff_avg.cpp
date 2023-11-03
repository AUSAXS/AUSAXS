#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <form_factor/FormFactor.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>
#include <table/DebyeTable.h>
#include <settings/All.h>
#include <constants/Constants.h>

using namespace data::record;
using namespace data;

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

// calculation: 8 points
//          1 line  of length 0
//          3 lines of length 2
//          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
//          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
//
// calculation: 1 center point
//          1 line  of length 0
//          16 lines of length sqrt(3) = 1.73 (counting both directions)
//
// sum:
//          9 line  of length 0
//          16 lines of length sqrt(3)
//          24 lines of length 2
//          24 lines of length sqrt(8)
//          8 lines of length sqrt(12)
auto width = constants::axes::d_axis.width();
std::vector<double> d = {
    0, 
    constants::axes::d_vals[std::round(std::sqrt(3)/width)], 
    constants::axes::d_vals[std::round(2./width)], 
    constants::axes::d_vals[std::round(std::sqrt(8)/width)], 
    constants::axes::d_vals[std::round(std::sqrt(12)/width)]
};

#define DEBYE_DEBUG 0
// unsigned int qcheck = 0;
TEST_CASE("CompositeDistanceHistogramFFAvg::debye_transform") {
    settings::molecule::use_effective_charge = false;
    auto ff_carbon = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::C);
    auto ff_exv = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::EXCLUDED_VOLUME);
    auto ff_w = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::OH);
    const auto& q_axis = constants::axes::q_vals;
    std::vector<double> Iq_exp(q_axis.size(), 0);

    SECTION("no water") {
        std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b5 = {Atom(Vector3<double>( 0,  0,  0), 1, constants::atom_t::C, "C", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
        Molecule protein(a);

        set_unity_charge(protein);
        double Z = protein.get_excluded_volume()*constants::charge::density::water/9;
        protein.set_excluded_volume_scaling(1./Z);

        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double aasum = 
                9 + 
                16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) +
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                8 *std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

            double axsum = 
                16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) +
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                8 *std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "aasum: " << aasum << std::endl;
                    std::cout << "\t9*1" << std::endl;
                    std::cout << "\t16*sin(" << q_axis[q] << "*" << d[1] << ")/(" << q_axis[q] << "*" << d[1] << ") = " << std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) << std::endl;
                    std::cout << "\t24*sin(" << q_axis[q] << "*" << d[2] << ")/(" << q_axis[q] << "*" << d[2] << ") = " << std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) << std::endl;
                    std::cout << "\t24*sin(" << q_axis[q] << "*" << d[3] << ")/(" << q_axis[q] << "*" << d[3] << ") = " << std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) << std::endl;
                    std::cout << "\t8*sin(" << q_axis[q] << "*" << d[4] << ")/(" << q_axis[q] << "*" << d[4] << ") = " << std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]) << std::endl;
                }
            #endif

            Iq_exp[q] += aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);
            Iq_exp[q] -= 2*axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);
            Iq_exp[q] += aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2);
        }

        auto Iq = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
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

        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double awsum = 8*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]);
            double aasum = 
                8 + 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
            double axsum = 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

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
            Iq_exp[q] -= 2*axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);  // -2ax
            Iq_exp[q] += aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2);                     // + xx
            Iq_exp[q] += 2*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);    // +2aw
            Iq_exp[q] -= 2*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);       // -2wx
            Iq_exp[q] += 1*std::pow(ff_w.evaluate(q_axis[q]), 2);                           // + ww

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "aasum = " << aasum << std::endl;
                    std::cout << "awsum = " << awsum << std::endl;
                    std::cout << "(aa) Iq_exp[" << q << "] += " << aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(ax) Iq_exp[" << q << "] -= " << 2*axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(aw) Iq_exp[" << q << "] += " << 2*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(xx) Iq_exp[" << q << "] += " << aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(wx) Iq_exp[" << q << "] -= " << 2*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(ww) Iq_exp[" << q << "] += " << 1*std::pow(ff_w.evaluate(q_axis[q]), 2) << std::endl;
                }
            #endif
        }
        auto Iq = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("real scalings") {
        std::vector<Atom> b1 =  {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 =  {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 =  {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 =  {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1),  Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Water> w = {Water(Vector3<double>( 0,  0,  0), 1, constants::atom_t::O, "HOH", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Molecule protein(a, w);

        double ZX = protein.get_excluded_volume()*constants::charge::density::water/8;
        double ZC = 6; // 6
        double ZO = 8; // 8

        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double awsum = 8*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]);
            double aasum = 
                8 + 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
            double axsum = 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "aasum = " << ZC*ZC*aasum << std::endl;
                    std::cout << "axsum = " << ZC*ZX*axsum << std::endl;
                    std::cout << "awsum = " << ZC*ZO*awsum << std::endl;
                    std::cout << "(aa) Iq_exp[" << q << "] += " << "Za*Za*aasum*ffa*ffa   =   " << ZC*ZC*aasum << "*" << ff_carbon.evaluate(q_axis[q])*ff_carbon.evaluate(q_axis[q]) << " = " << ZC*ZC*aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(ax) Iq_exp[" << q << "] -= " << "2*Zx*Zx*axsum*ffa*ffx = 2*" << ZX*ZX*axsum << "*" << ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << " = " << 2*ZC*ZX*axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(aw) Iq_exp[" << q << "] += " << "2*Zx*Zw*awsum*ffa*ffw = 2*" << ZX*ZO*awsum << "*" << ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << " = " << 2*ZX*ZO*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(xx) Iq_exp[" << q << "] += " << "Zx*Zx*xxsum*ffx*ffx   =   " << ZX*ZX*aasum << "*" << ff_exv.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << " = "<< ZX*ZX*aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(wx) Iq_exp[" << q << "] -= " << "2*Zx*Zw*wxsum*ffx*ffw = 2*" << ZX*ZO*awsum << "*" << ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << " = " << 2*ZX*ZO*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(ww) Iq_exp[" << q << "] += " << "Zw*Zw*wwsum*ffw*ffw   =   " << ZO*ZO*1 << "*" << ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << " = " << 1*ZO*ZO*std::pow(ff_w.evaluate(q_axis[q]), 2) << std::endl;
                }
            #endif

            Iq_exp[q] += ZC*ZC*aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);                    // + aa
            Iq_exp[q] -= 2*ZC*ZX*axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);    // -2ax
            Iq_exp[q] += ZX*ZX*aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2);                       // + xx
            Iq_exp[q] += 2*ZC*ZO*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);      // +2aw
            Iq_exp[q] -= 2*ZO*ZX*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);         // -2wx
            Iq_exp[q] += 1*ZO*ZO*std::pow(ff_w.evaluate(q_axis[q]), 2);                             // + ww
        }
        auto Iq = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("real data") {
        Molecule protein("test/files/2epe.pdb");
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
        auto Iq = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();

        unsigned int N = protein.atom_size();
        unsigned int M = protein.water_size();
        REQUIRE_THAT(Iq[0], Catch::Matchers::WithinRel(std::pow(N*ZC + M*ZO - N*ZX, 2) + 2*N*ZC*ZX, 1e-3));
    }
}

TEST_CASE("CompositeDistanceHistogramFFAvg::get_profile") {
    data::Molecule protein("test/files/2epe.pdb");
    auto hist_data = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all();
    auto hist = static_cast<hist::CompositeDistanceHistogramFFAvg*>(hist_data.get());
    auto Iq = hist->debye_transform();
    auto profile_sum = 
          hist->get_profile_ww() - hist->get_profile_wx() + hist->get_profile_xx()
        + hist->get_profile_aw() - hist->get_profile_ax() + hist->get_profile_aa();
    REQUIRE(Iq.size() == profile_sum.size());
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        REQUIRE_THAT(Iq[i], Catch::Matchers::WithinAbs(profile_sum[i], 1e-3));
    }
}