#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <form_factor/NormalizedFormFactor.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>
#include <table/DebyeTable.h>
#include <settings/All.h>
#include <constants/Constants.h>

#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace data;

#define DEBYE_DEBUG 0
// unsigned int qcheck = 0;
TEST_CASE("CompositeDistanceHistogramFFAvg::debye_transform") {
    settings::molecule::implicit_hydrogens = false;
    auto ff_carbon = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C);
    auto ff_exv = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::EXCLUDED_VOLUME);
    auto ff_w = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::OH);
    const auto& q_axis = constants::axes::q_vals;
    std::vector<double> Iq_exp(q_axis.size(), 0);
    auto d = SimpleCube::d;

    SECTION("no water") {
        std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b5 = {AtomFF({ 0,  0,  0}, form_factor::form_factor_t::C)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
        DebugMolecule protein(a);

        set_unity_charge(protein);
        double Z = protein.get_volume_grid()*constants::charge::density::water/9;
        protein.set_volume_scaling(1./Z);

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

                    std::cout << "axsum: " << axsum << std::endl;
                    std::cout << "\t16*sin(" << q_axis[q] << "*" << d[1] << ")/(" << q_axis[q] << "*" << d[1] << ") = " << std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) << std::endl;
                    std::cout << "\t24*sin(" << q_axis[q] << "*" << d[2] << ")/(" << q_axis[q] << "*" << d[2] << ") = " << std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) << std::endl;
                    std::cout << "\t24*sin(" << q_axis[q] << "*" << d[3] << ")/(" << q_axis[q] << "*" << d[3] << ") = " << std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) << std::endl;
                    std::cout << "\t8*sin(" << q_axis[q] << "*" << d[4] << ")/(" << q_axis[q] << "*" << d[4] << ") = " << std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]) << std::endl;

                    std::cout << "Iaa: " << aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "Iax: " << axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]) << std::endl;
                    std::cout << "Ixx: " << aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2) << std::endl;
                }
            #endif

            Iq_exp[q] += aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);
            Iq_exp[q] -= 2*axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);
            Iq_exp[q] += aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2);
        }

        auto Iq = hist::HistogramManagerMTFFAvg<false, false>(&protein).calculate_all()->debye_transform();
        CHECK(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("with water") {
        std::vector<AtomFF> b1 =  {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C),  AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 =  {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C),  AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 =  {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C),  AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 =  {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C),  AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<Water>  w  =  {Water({0,  0,  0})};
        std::vector<Body> a = {Body(b1, w), Body(b2), Body(b3), Body(b4)};
        DebugMolecule protein(a);

        set_unity_charge(protein);
        double Z = protein.get_volume_grid()*constants::charge::density::water/8;
        protein.set_volume_scaling(1./Z);

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
        auto Iq = hist::HistogramManagerMTFFAvg<false, false>(&protein).calculate_all()->debye_transform();
        CHECK(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("real scalings") {
        std::vector<AtomFF> b1 =  {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 =  {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 =  {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 =  {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<Water>  w  =  {Water({0,  0,  0})};
        std::vector<Body> a = {Body(b1, w), Body(b2), Body(b3), Body(b4)};
        DebugMolecule protein(a);

        double ZX = protein.get_volume_grid()*constants::charge::density::water/8;

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

            Iq_exp[q] += aasum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);                    // + aa
            Iq_exp[q] -= 2*ZX*axsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]); // -2ax
            Iq_exp[q] += ZX*ZX*aasum*std::pow(ff_exv.evaluate(q_axis[q]), 2);                 // + xx
            Iq_exp[q] += 2*awsum*ff_carbon.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);      // +2aw
            Iq_exp[q] -= 2*awsum*ff_exv.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);         // -2wx
            Iq_exp[q] += 1*std::pow(ff_w.evaluate(q_axis[q]), 2);                             // + ww
        }
        auto Iq = hist::HistogramManagerMTFFAvg<false, false>(&protein).calculate_all()->debye_transform();
        CHECK(compare_hist(Iq_exp, Iq.get_counts()));
    }

    // TODO: fix this test
    // SECTION("real data") {
    //     settings::general::verbose = false;
    //     DebugMolecule protein("tests/files/2epe.pdb");
    //     double ZX = protein.get_volume_grid()*constants::charge::density::water/protein.atom_size();
    //     double ZC = 6;
    //     double ZO = 8;
    //     for (auto& body : protein.get_bodies()) {
    //         for (auto& atom : body.get_atoms()) {
    //             atom.set_element(constants::atom_t::C);
    //             atom.set_effective_charge(ZC);
    //         }
    //     }
    //     protein.generate_new_hydration();
    //     for (auto& water : protein.get_waters()) {
    //         water.set_effective_charge(ZO);
    //     }
    //     auto Iq = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();

    //     unsigned int N = protein.atom_size();
    //     unsigned int M = protein.water_size();
    //     REQUIRE_THAT(Iq[0], Catch::Matchers::WithinRel(std::pow(N*ZC + M*ZO - N*ZX, 2) + 2*N*ZC*ZX, 1e-3));
    // }
}

TEST_CASE("CompositeDistanceHistogramFFAvg::get_profile") {
    settings::general::verbose = false;
    data::Molecule protein("tests/files/2epe.pdb");
    auto hist_data = hist::HistogramManagerMTFFAvg<false, false>(&protein).calculate_all();
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

TEST_CASE("CompositeDistanceHistogramFFAvg: Debye-Waller factors") {
    settings::general::verbose = false;
    data::Molecule protein("tests/files/2epe.pdb");
    auto h = hist::HistogramManagerMTFFAvg<false, false>(&protein).calculate_all();
    auto h_cast = static_cast<hist::CompositeDistanceHistogramFFAvg*>(h.get());

    auto analytical = [&] (double Ba_, double Bx_) {
        h_cast->apply_atomic_debye_waller_factor(0);
        h_cast->apply_exv_debye_waller_factor(0);
        auto aa = h_cast->get_profile_aa();
        auto aw = h_cast->get_profile_aw();
        auto ww = h_cast->get_profile_ww();
        auto xx = h_cast->get_profile_xx();
        auto ax = h_cast->get_profile_ax();
        auto wx = h_cast->get_profile_wx();
        auto Iq = std::vector<double>(aa.size(), 0);
        for (unsigned int i = 0; i < Iq.size(); ++i) {
            double Ba = hist::CompositeDistanceHistogramFFAvg::get_atomic_debye_waller_factor(constants::axes::q_vals[i], Ba_);
            double Bx = hist::CompositeDistanceHistogramFFAvg::get_exv_debye_waller_factor(constants::axes::q_vals[i], Bx_);
            Iq[i] = Ba*Ba*aa[i] + Ba*aw[i] + ww[i] + Bx*Bx*xx[i] - Ba*Bx*ax[i] - Bx*wx[i];
        }
        return Iq;
    };

    SECTION("no change") {
        auto Iq_base = analytical(0, 0);
        h_cast->apply_atomic_debye_waller_factor(0);
        h_cast->apply_exv_debye_waller_factor(0);
        auto Iq = h_cast->debye_transform();
        REQUIRE(compare_hist(Iq_base, Iq));
    }

    SECTION("atomic") {
        auto B = GENERATE(0.1, 0.5, 1.0, 2.0);
        auto Iq_exp = analytical(B, 0);
        h_cast->apply_atomic_debye_waller_factor(B);
        auto Iq = h_cast->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq));
    }

    SECTION("excluded volume") {
        auto B = GENERATE(0.1, 0.5, 1.0, 2.0);
        auto Iq_exp = analytical(0, B);
        h_cast->apply_exv_debye_waller_factor(B);
        auto Iq = h_cast->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq));
    }

    SECTION("both") {
        auto Ba = GENERATE(0.1, 0.5, 1.0, 2.0);
        auto Bx = GENERATE(0.1, 0.5, 1.0, 2.0);
        auto Iq_exp = analytical(Ba, Bx);
        h_cast->apply_atomic_debye_waller_factor(Ba);
        h_cast->apply_exv_debye_waller_factor(Bx);
        auto Iq = h_cast->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq));
    }
}