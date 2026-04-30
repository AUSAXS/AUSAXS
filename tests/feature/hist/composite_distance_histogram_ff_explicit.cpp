#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>
#include <table/DebyeTable.h>
#include <settings/All.h>
#include <constants/Constants.h>

#include "hist/hist_test_helper.h"
#include "settings/GeneralSettings.h"

using namespace ausaxs;
using namespace ausaxs::data;

#define DEBYE_DEBUG 0
unsigned int qcheck = 26;
TEST_CASE("CompositeDistanceHistogramFFAvg::debye_transform") {
    settings::molecule::implicit_hydrogens = false;
    auto ff_C = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C);
    auto ff_w = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::OH);
    auto ff_Cx = form_factor::ExvTableManager::get_current_exv_form_factor_set().get(form_factor::form_factor_t::C);
    auto ff_wx = form_factor::ExvTableManager::get_current_exv_form_factor_set().get(form_factor::form_factor_t::OH);
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

            Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);
            Iq_exp[q] -= 2*aasum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);
            Iq_exp[q] += aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2);
        }

        auto Iq = hist::HistogramManagerMTFFExplicit<false, false>(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("with water") {
        std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<Water> w = {Water({0,  0,  0})};
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

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "8*1" << std::endl;
                    std::cout << "8*sin(" << q_axis[q] << "*" << d[4] << ")/(" << q_axis[q] << "*" << d[4] << ") = " << std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]) << std::endl;
                    std::cout << "24*sin(" << q_axis[q] << "*" << d[1] << ")/(" << q_axis[q] << "*" << d[1] << ") = " << std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) << std::endl;
                    std::cout << "24*sin(" << q_axis[q] << "*" << d[2] << ")/(" << q_axis[q] << "*" << d[2] << ") = " << std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) << std::endl;
                    std::cout << "8*sin(" << q_axis[q] << "*" << d[3] << ")/(" << q_axis[q] << "*" << d[3] << ") = " << std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) << std::endl;
                }
            #endif
            Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);                 // + aa
            Iq_exp[q] -= 2*aasum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);  // -2ax
            Iq_exp[q] += aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2);                // + xx
            Iq_exp[q] += 2*awsum*ff_C.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);   // +2aw
            Iq_exp[q] -= awsum*(
                            ff_w.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]) 
                            + ff_wx.evaluate(q_axis[q])*ff_C.evaluate(q_axis[q])
                        );  // -wx-xw
            Iq_exp[q] += 1*std::pow(ff_w.evaluate(q_axis[q]), 2);                     // + ww

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "aasum = " << aasum << std::endl;
                    std::cout << "awsum = " << awsum << std::endl;
                    std::cout << "(aa) Iq_exp[" << q << "] += " << aasum*std::pow(ff_C.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(ax) Iq_exp[" << q << "] -= " << 2*axsum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(xx) Iq_exp[" << q << "] += " << aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(aw) Iq_exp[" << q << "] += " << 2*awsum*ff_C.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(wx) Iq_exp[" << q << "] -= " << 2*awsum*ff_w.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(ww) Iq_exp[" << q << "] += " << 1*std::pow(ff_w.evaluate(q_axis[q]), 2) << std::endl;
                }
            #endif
        }
        auto Iq = hist::HistogramManagerMTFFExplicit<false, false>(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    SECTION("real scalings") {
        std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<Water> w = {Water({0,  0,  0})};
        std::vector<Body> a = {Body(b1, w), Body(b2), Body(b3), Body(b4)};
        Molecule protein(a);

        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double awsum = 8*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]);
            double aasum = 
                8 + 
                24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

            #if DEBYE_DEBUG
                if (q==qcheck) {
                    std::cout << "aasum = " << aasum << std::endl;
                    std::cout << "axsum = " << axsum << std::endl;
                    std::cout << "awsum = " << awsum << std::endl;
                    std::cout << "(aa) Iq_exp[" << q << "] += " << ZC*ZC*aasum*std::pow(ff_C.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(ax) Iq_exp[" << q << "] -= " << 2*ZC*ZX*axsum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(aw) Iq_exp[" << q << "] += " << 2*ZX*ZO*awsum*ff_C.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(xx) Iq_exp[" << q << "] += " << ZX*ZX*aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2) << std::endl;
                    std::cout << "(xw) Iq_exp[" << q << "] -= " << 2*ZX*ZO*awsum*ff_wx.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]) << std::endl;
                    std::cout << "(ww) Iq_exp[" << q << "] += " << 1*ZO*ZO*std::pow(ff_w.evaluate(q_axis[q]), 2) << std::endl;
                }
            #endif

            Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);                 // + aa
            Iq_exp[q] -= 2*aasum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);  // -2ax
            Iq_exp[q] += aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2);                // + xx
            Iq_exp[q] += 2*awsum*ff_C.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);   // +2aw
            Iq_exp[q] -= awsum*(
                            ff_w.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q])
                            + ff_wx.evaluate(q_axis[q])*ff_C.evaluate(q_axis[q])
                        );  // -2wx
            Iq_exp[q] += 1*std::pow(ff_w.evaluate(q_axis[q]), 2);                     // + ww
        }
        auto Iq = hist::HistogramManagerMTFFExplicit<false, false>(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    // SECTION("real data") {
    //     Molecule protein("tests/files/2epe.pdb");
    //     double ZX = protein.get_excluded_volume()*constants::charge::density::water/protein.atom_size();
    //     double ZC = 6;
    //     double ZO = 8;
    //     for (auto& body : protein.get_bodies()) {
    //         for (auto& atom : body.get_atoms()) {
    //             atom.set_element(form_factor::form_factor_t::C);
    //             atom.set_effective_charge(ZC);
    //         }
    //     }
    //     protein.generate_new_hydration();
    //     for (auto& water : protein.get_waters()) {
    //         water.set_effective_charge(ZO);
    //     }
    //     auto Iq = hist::HistogramManagerMTFFExplicit(&protein).calculate_all()->debye_transform();

    //     unsigned int N = protein.atom_size();
    //     unsigned int M = protein.water_size();
    //     REQUIRE_THAT(Iq[0], Catch::Matchers::WithinRel(std::pow(N*ZC + M*ZO - N*ZX, 2) + 2*N*ZC*ZX, 1e-3));
    // }
}

TEST_CASE("CompositeDistanceHistogramFFAvg::get_profile") {
    settings::general::verbose = false;

    data::Molecule protein("tests/files/2epe.pdb");
    auto hist_data = hist::HistogramManagerMTFFExplicit<false, false>(&protein).calculate_all();
    auto hist = static_cast<hist::CompositeDistanceHistogramFFExplicit*>(hist_data.get());
    auto Iq = hist->debye_transform();
    auto profile_sum = 
          hist->get_profile_ww() - hist->get_profile_wx() + hist->get_profile_xx()
        + hist->get_profile_aw() - hist->get_profile_ax() + hist->get_profile_aa();
    REQUIRE(Iq.size() == profile_sum.size());
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        REQUIRE_THAT(Iq[i], Catch::Matchers::WithinRel(profile_sum[i], 1e-3));
    }
}

// #include <dataset/SimpleDataset.h>
// #include <plots/PlotDataset.h>
// struct DummyCDHFFX : public hist::CompositeDistanceHistogramFFExplicit {
//     double Gq(double q) const {return exv_factor(q);}
// };
// TEST_CASE("plot_Gq", "[manual]") {
//     SimpleDataset Gq;
//     double r0 = 1.5;
//     DummyCDHFFX hist;
//     hist.apply_excluded_volume_scaling_factor(r0);
//     for (unsigned int i = 0; i < constants::axes::q_vals.size(); ++i) {
//         Gq.push_back(constants::axes::q_vals[i], hist.Gq(constants::axes::q_vals[i]));
//     }
//     plots::PlotDataset::quick_plot(Gq, plots::PlotOptions(), io::File("temp/tests/composite_distance_histogram_ff_explicit/Gq.png"));
// }