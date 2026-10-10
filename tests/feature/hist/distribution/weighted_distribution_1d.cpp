#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <settings/ExvSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/MoleculeSettings.h>

#include <support/hist_test_helper.h>

#include <algorithm>
#include <numeric>
#include <ranges>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::data;

class DistanceHistogramDebug : public DistanceHistogram {
    public:
        DistanceHistogramDebug(DistanceHistogram&& other) : DistanceHistogram(std::move(other)) {}
        auto get_sinc_table() const {return sinc_table.get_sinc_table();}
};

TEST_CASE("WeightedDistribution: sinc_table") {
    settings::molecule::implicit_hydrogens = false;
    std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b5 = {AtomFF({ 0,  0,  0}, form_factor::form_factor_t::C)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
    Molecule protein(a);

    auto hist = hist::HistogramManagerMT<true>(&protein).calculate_all();
    auto Iq = hist->debye_transform();

    const auto& bins = constants::axes::d_vals;
    const auto* table = DistanceHistogramDebug(std::move(hist)).get_sinc_table();
    for (int q = 0; q < static_cast<int>(table->size_q()); ++q) {
        std::vector<double> sinc(20);
        for (int d = 0; d < 20; ++d) {
            double qd = constants::axes::q_vals[q]*bins[d];
            double val = 0;
            if (qd < 1e-3) {val = 1 - qd*qd/6 + qd*qd*qd*qd/120;}
            else {val = std::sin(qd)/qd;}
            REQUIRE_THAT(table->lookup(q, d), Catch::Matchers::WithinAbs(val, 1e-6));
            sinc[d] = val;
        }
        std::ranges::transform(sinc, std::ranges::subrange(table->begin(q), table->end(q)), sinc.begin(), std::minus<>());
        REQUIRE_THAT(std::reduce(sinc.begin(), sinc.end(), 0.0), Catch::Matchers::WithinAbs(0, 1e-6));
    }
}

// Check that the weighted distance axis is correctly calculated for all histogram managers.
TEST_CASE("WeightedDistribution: distance_calculators") {
    settings::molecule::implicit_hydrogens = false;
    std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b5 = {AtomFF({ 0,  0,  0}, form_factor::form_factor_t::C)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
    Molecule protein(a);
    form_factor::manager::use_form_factors(protein); // the form factor managers below are constructed directly, bypassing the factory which normally selects these

    { // hm
        CHECK(SimpleCube::check_default(hist::HistogramManager<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManager<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt
        CHECK(SimpleCube::check_default(hist::HistogramManagerMT<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMT<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt_ff_avg
        CHECK(SimpleCube::check_default(hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMTFFAvg<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt_ff_explicit
        CHECK(SimpleCube::check_default(hist::HistogramManagerMTFFExplicit<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMTFFExplicit<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // hm_mt_ff_grid
        CHECK(SimpleCube::check_exact(hist::HistogramManagerMTFFGrid(&protein).calculate_all()->get_d_axis()));
        // CHECK(check_exact(hist::HistogramManagerMTFFGrid(&protein).calculate_all())); // exv cells dominates bin locs in this case
    }
    { // phm
        CHECK(SimpleCube::check_default(hist::PartialHistogramManager<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::PartialHistogramManager<true>(&protein).calculate_all()->get_d_axis()));
    }
    { // phm_mt
        CHECK(SimpleCube::check_default(hist::PartialHistogramManagerMT<false>(&protein).calculate_all()->get_d_axis()));
        CHECK(SimpleCube::check_exact(hist::PartialHistogramManagerMT<true>(&protein).calculate_all()->get_d_axis()));
    }
}

// Check that the basic histogram managers agree on a weighted debye transform.
TEST_CASE("CompositeDistanceHistogram::debye_transform (weighted)") {
    settings::exv::exv_method = settings::exv::ExvMethod::None; // the expected histograms use the unmodified atomic weights
    settings::molecule::implicit_hydrogens = false;
    settings::general::warnings = true;
    auto d_exact = SimpleCube::d_exact;

    SECTION("no water") {
        std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b5 = {AtomFF({ 0,  0,  0}, form_factor::form_factor_t::C)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4), Body(b5)};
        Molecule protein(a);

        set_unity_charge(protein);

        std::vector<double> Iq_exp;
        {
            const auto& q_axis = constants::axes::q_vals;
            Iq_exp.resize(q_axis.size(), 0);
            auto ff2 = [] (double q) {return std::exp(-q*q);};

            for (int q = 0; q < static_cast<int>(q_axis.size()); ++q) {
                double dsum = 
                    9 + 
                    16*std::sin(q_axis[q]*d_exact[1])/(q_axis[q]*d_exact[1]) +
                    24*std::sin(q_axis[q]*d_exact[2])/(q_axis[q]*d_exact[2]) + 
                    24*std::sin(q_axis[q]*d_exact[3])/(q_axis[q]*d_exact[3]) + 
                    8 *std::sin(q_axis[q]*d_exact[4])/(q_axis[q]*d_exact[4]);
                Iq_exp[q] += dsum*ff2(q_axis[q]);
            }
        }

        {
            auto Iq = hist::HistogramManager<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_intensity()));
        }
        {
            auto Iq = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_intensity()));
        }
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

        std::vector<double> Iq_exp;
        {
            const auto& q_axis = constants::axes::q_vals;
            Iq_exp.resize(q_axis.size(), 0);
            auto ff = [] (double q) {return std::exp(-q*q/2);};

            for (int q = 0; q < static_cast<int>(q_axis.size()); ++q) {
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
            auto Iq = hist::HistogramManager<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_intensity()));
        }
        {
            auto Iq = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();
            REQUIRE(compare_hist(Iq_exp, Iq.get_intensity()));
        }
    }
}
