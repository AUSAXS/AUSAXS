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

using namespace ausaxs;
using namespace ausaxs::data;

auto test = [] (data::Molecule& protein, auto&& phm) {
    // no changes
    auto p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->debye_transform();
    auto phm_res = phm(protein)->debye_transform();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // add symmetry
    // protein.get_body(0).symmetry().add({Vector3<double>(1, 0, 0)});
    // p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->debye_transform();
    // phm_res = phm(protein)->debye_transform();
    // REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry
    protein.get_body(0).symmetry().get(0).translate = {0, 1, 0};
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->debye_transform();
    phm_res = phm(protein)->debye_transform();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));
};

// Test that subsequent calculations are correct
TEST_CASE("PartialSymmetryManagerMT: subsequent calculations") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    data::Molecule protein({
        Body("tests/files/2epe.pdb"), 
        Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
    });

    protein.get_body(0).symmetry().add({Vector3<double>(1, 0, 0)});
    protein.generate_new_hydration();
    test(protein, [](const Molecule& protein) {return hist::PartialSymmetryManagerMT<true>(&protein).calculate_all();});
}