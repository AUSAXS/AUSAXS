#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/Symmetry.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "settings/GeneralSettings.h"
#include "settings/HistogramSettings.h"

using namespace ausaxs;
using namespace ausaxs::data;

auto test = [] (data::Molecule& protein) {
    // no changes
    auto p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    auto phm_res = protein.get_histogram()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // add symmetry
    // protein.get_body(0).symmetry().add({Vector3<double>(1, 0, 0)});
    // p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->debye_transform();
    // phm_res = phm(protein)->debye_transform();
    // REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry
    protein.get_body(0).symmetry().get(0).translate = {0, 1, 0};
    phm_res = protein.get_histogram()->get_total_counts();
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));
};

// Test that subsequent calculations are correct
TEST_CASE("PartialSymmetryManagerMT: subsequent calculations") {
    settings::general::threads = 1;
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT;

    SECTION("simple") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}}
        });
        set_unity_charge(protein);
        protein.get_body(0).symmetry().add({{-1, 0, 0}});
        test(protein);
    }

    // SECTION("2epe") {
    //     data::Molecule protein({
    //         Body("tests/files/2epe.pdb"), 
    //         Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
    //     });
    //     protein.get_body(0).symmetry().add({Vector3<double>(1, 0, 0)});
    //     protein.generate_new_hydration();    
    //     test(protein);
    // }
}