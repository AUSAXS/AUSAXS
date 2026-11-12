#include <algorithm>
#include <utility>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset.h>
#include <form_factor/xray/XrayFormFactor.h>
#include <hist/scattering_model/ScatteringModelPoint.h>
#include <settings/All.h>

#include <support/hist_test_helper.h>


using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("Wide-angle scattering: gold nanoparticles") {
    settings::general::input_q_unit = settings::general::QUnit::USER_A;
    settings::axes::clamp_to_qrange = false;
    settings::axes::bin_width = 0.1;

    // expected scattering
    Dataset expected("tests/files/carbon_sphere_waxs.dat");

    // calculate scattering from the molecule
    Molecule mol("tests/files/carbon_sphere.xyz");
    mol.clear_hydration();

    auto normalise = [] (std::vector<double> y) {
        double s = y.front();
        std::ranges::transform(y, y.begin(), [s] (double v) {return v/s;});
        return y;
    };

    SECTION("point model") {
        settings::exv::exv_method = settings::exv::ExvMethod::None;
        auto I = hist::ScatteringModelPoint(std::move(*mol.get_total_histogram())).debye_transform(expected.x());
        REQUIRE(compare_hist_approx(normalise(I.y()), normalise(expected.y()), 1e-6, 0.01));
    }

    // TODO: a placeholder until the exv models can be evaluated at arbitrary q, at which point they can be compared with the vacuum model there.
    // every atom is a carbon, so the form factor is an overall factor f_C(q)^2 which is divided out
    SECTION("vacuum model") {
        REQUIRE(std::ranges::all_of(mol.get_atoms(), [] (const AtomFF& a) {return a.form_factor_type() == form_factor::form_factor_t::C;}));
        settings::exv::exv_method = settings::exv::ExvMethod::Vacuum;
        mol.reset_histogram_manager();
        auto I = mol.get_histogram()->debye_transform(expected.x());
        auto ff = form_factor::xray::raw::get(form_factor::form_factor_t::C);
        std::ranges::transform(I.y(), I.x(), I.y().begin(), [&ff] (double Iq, double q) {return Iq/ff.evaluate_self(q);});
        REQUIRE(compare_hist_approx(normalise(I.y()), normalise(expected.y()), 1e-6, 0.01));
    }
}
