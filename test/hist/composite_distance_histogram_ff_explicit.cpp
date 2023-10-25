#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/CompositeDistanceHistogramFFAvg.h>
#include <hist/CompositeDistanceHistogramFFExplicit.h>
#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
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

#define DEBYE_DEBUG 0
unsigned int qcheck = 26;
TEST_CASE("CompositeDistanceHistogramFFAvg::debye_transform") {
    settings::molecule::use_effective_charge = false;
    auto ff_C = form_factor::storage::get_form_factor(form_factor::form_factor_t::C);
    auto ff_w = form_factor::storage::get_form_factor(form_factor::form_factor_t::O);
    auto ff_Cx = form_factor::storage::exv::get_form_factor(form_factor::form_factor_t::C);
    // auto ff_wx = form_factor::storage::exv::get_form_factor(form_factor::form_factor_t::O);
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

            Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);
            Iq_exp[q] -= 2*axsum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);
            Iq_exp[q] += aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2);
        }

        auto Iq = hist::HistogramManagerMTFFExplicit(&protein).calculate_all()->debye_transform();
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
            Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);                 // + aa
            Iq_exp[q] -= 2*axsum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);  // -2ax
            Iq_exp[q] += aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2);                // + xx
            Iq_exp[q] += 2*awsum*ff_C.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);   // +2aw
            Iq_exp[q] -= 2*awsum*ff_w.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);  // -2wx
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
        auto Iq = hist::HistogramManagerMTFFExplicit(&protein).calculate_all()->debye_transform();
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

            Iq_exp[q] += ZC*ZC*aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);                 // + aa
            Iq_exp[q] -= 2*ZC*axsum*ff_C.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);     // -2ax
            Iq_exp[q] += aasum*std::pow(ff_Cx.evaluate(q_axis[q]), 2);                      // + xx
            Iq_exp[q] += 2*ZC*ZO*awsum*ff_C.evaluate(q_axis[q])*ff_w.evaluate(q_axis[q]);   // +2aw
            Iq_exp[q] -= 2*ZO*awsum*ff_w.evaluate(q_axis[q])*ff_Cx.evaluate(q_axis[q]);     // -2wx
            Iq_exp[q] += 1*ZO*ZO*std::pow(ff_w.evaluate(q_axis[q]), 2);                     // + ww
        }
        auto Iq = hist::HistogramManagerMTFFExplicit(&protein).calculate_all()->debye_transform();
        REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
    }

    // SECTION("real data") {
    //     Molecule protein("test/files/2epe.pdb");
    //     double ZX = protein.get_excluded_volume()*constants::charge::density::water/protein.atom_size();
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
    //     auto Iq = hist::HistogramManagerMTFFExplicit(&protein).calculate_all()->debye_transform();

    //     unsigned int N = protein.atom_size();
    //     unsigned int M = protein.water_size();
    //     REQUIRE_THAT(Iq[0], Catch::Matchers::WithinRel(std::pow(N*ZC + M*ZO - N*ZX, 2) + 2*N*ZC*ZX, 1e-3));
    // }
}

TEST_CASE("CompositeDistanceHistogramFFAvg::get_profile") {
    data::Molecule protein("test/files/2epe.pdb");
    auto hist_data = hist::HistogramManagerMTFFExplicit(&protein).calculate_all();
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

#include <dataset/SimpleDataset.h>
#include <plots/PlotDataset.h>
struct DummyCDHFFX : public hist::CompositeDistanceHistogramFFExplicit {
    double Gq(double q) const {return G_factor(q);}
};
TEST_CASE("plot_Gq", "[manual]") {
    SimpleDataset Gq;
    double r0 = 1.5;
    DummyCDHFFX hist;
    hist.apply_excluded_volume_scaling_factor(r0);
    for (unsigned int i = 0; i < constants::axes::q_vals.size(); ++i) {
        Gq.push_back(constants::axes::q_vals[i], hist.Gq(constants::axes::q_vals[i]));
    }
    plots::PlotDataset::quick_plot(Gq, io::File("temp/test/composite_distance_histogram_ff_explicit/Gq.png"));
}

TEST_CASE("plot_cmp", "[manual]") {
    data::Molecule protein("data/SASDJG5/SASDJG5.pdb");
    protein.generate_new_hydration();
    auto h = hist::HistogramManagerMTFFAvg(&protein).calculate_all();
    auto hx = hist::HistogramManagerMTFFExplicit(&protein).calculate_all();
    auto hc = dynamic_cast<hist::CompositeDistanceHistogramFFAvg*>(h.get());
    auto hcx = dynamic_cast<hist::CompositeDistanceHistogramFFExplicit*>(hx.get());

    auto aa = hc->get_profile_aa().as_dataset();
    auto ax = hc->get_profile_ax().as_dataset();
    auto xx = hc->get_profile_xx().as_dataset();
    auto aw = hc->get_profile_aw().as_dataset();
    auto wx = hc->get_profile_wx().as_dataset();
    auto ww = hc->get_profile_ww().as_dataset();

    auto aa_x = hcx->get_profile_aa().as_dataset();
    auto ax_x = hcx->get_profile_ax().as_dataset();
    auto xx_x = hcx->get_profile_xx().as_dataset();
    auto aw_x = hcx->get_profile_aw().as_dataset();
    auto wx_x = hcx->get_profile_wx().as_dataset();
    auto ww_x = hcx->get_profile_ww().as_dataset();

    for (unsigned int i = 0; i < aa.size(); ++i) {
        aa.y(i) /= aa_x.y(i);
        ax.y(i) /= ax_x.y(i);
        xx.y(i) /= xx_x.y(i);
        aw.y(i) /= aw_x.y(i);
        wx.y(i) /= wx_x.y(i);
        ww.y(i) /= ww_x.y(i);
    }

    // double c_aa = aa.normalize();
    // double c_ax = ax.normalize();
    // double c_xx = xx.normalize();
    // double c_aw = aw.normalize();
    // double c_wx = wx.normalize();
    // double c_ww = ww.normalize();

    // convert number to scientific notation
    // std::stringstream ss_aa, ss_ax, ss_xx, ss_aw, ss_wx, ss_ww;
    // ss_aa << std::scientific << c_aa;
    // ss_ax << std::scientific << c_ax;
    // ss_xx << std::scientific << c_xx;
    // ss_aw << std::scientific << c_aw;
    // ss_wx << std::scientific << c_wx;
    // ss_ww << std::scientific << c_ww;

    // aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$c_{aa} = $" + ss_aa.str()}, {"color", style::color::red}});
    // ax.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$c_{ax} = $" + ss_ax.str()}, {"color", style::color::blue}});
    // xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$c_{xx} = $" + ss_xx.str()}, {"color", style::color::green}});
    // aw.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$c_{aw} = $" + ss_aw.str()}, {"color", style::color::pink}});
    // wx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$c_{wx} = $" + ss_wx.str()}, {"color", style::color::purple}});
    // ww.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$c_{ww} = $" + ss_ww.str()}, {"color", style::color::brown}});

    aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$I_{aa}$"}, {"color", style::color::red}});
    ax.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$I_{ax}$"}, {"color", style::color::blue}});
    xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$I_{xx}$"}, {"color", style::color::green}});
    aw.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$I_{aw}$"}, {"color", style::color::pink}});
    wx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$I_{wx}$"}, {"color", style::color::purple}});
    ww.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", false}, {"logy", true}, {"ylabel", "I"}, {"legend", "$I_{ww}$"}, {"color", style::color::brown}});
    plots::PlotDataset(aa)
        .plot(ax)
        .plot(xx)
        .plot(aw)
        .plot(wx)
        .plot(ww)
    .save("temp/test/composite_distance_histogram_ff_explicit/compare_ausaxs_x.png");
}