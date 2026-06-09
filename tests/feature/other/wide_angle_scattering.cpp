#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/SimpleDataset.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <plots/PlotDataset.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("Wide-angle scattering: gold nanoparticles") {
    settings::general::input_q_unit = settings::general::QUnit::USER_A;
    settings::axes::clamp_to_qrange = false;
    settings::exv::exv_method = settings::exv::ExvMethod::None;
    settings::axes::bin_width = 0.1;

    // expected scattering
    Dataset expected("tests/files/SphericalAuScattering.dat");

    // calculate scattering from the molecule
    Molecule mol("tests/files/SphericalAuNP.xyz");
    mol.clear_hydration();
    auto I = mol.get_total_histogram()->debye_transform(expected.x());
    std::transform(I.y().begin(), I.y().end(), I.x().begin(), I.y().begin(), [] (double I, double q) {
        return I*std::exp(q*q); // remove form factor imposed by debye transform
    });

    // compare the results
    plots::PlotDataset()
        .plot(expected, {{"legend", "expected"}, {"xlabel", "q"}, {"ylabel", "I(q)"}, {"logx", true}, {"logy", true}, {"color", "tab:orange"}})
        .plot(I, {{"legend", "calculated"}, {"color", "tab:blue"}})
    .save("temp/WideAngleScattering.png");
    REQUIRE(compare_hist_approx(I.y(), expected.y(), 1e-6, 0.01));
}