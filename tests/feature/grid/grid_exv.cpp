#include <catch2/catch_test_macros.hpp>

#include <grid/Grid.h>
#include <grid/exv/RawGridExv.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;

// check that the RawGridExv and RawGridWithSurfaceExv are consistent
TEST_CASE("RawGridExv: consistency") {
    settings::general::verbose = false;
    settings::grid::min_bins = 100;

    SECTION("simple") {
        std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Molecule protein(a);
        auto grid = protein.get_grid();

        auto vol = grid::exv::RawGridExv::create(grid);
        auto vol2 = grid::exv::RawGridWithSurfaceExv::create(grid, false);

        REQUIRE(vol.interior.size() == vol2.interior.size());
        REQUIRE(vol.surface.size() == vol2.surface.size());

        for (unsigned int i = 0; i < vol.interior.size(); i++) {
            CHECK(vol.interior[i] == vol2.interior[i]);
        }
    }

    SECTION("real data") {
        settings::general::verbose = false;
        Molecule protein("tests/files/2epe.pdb");
        protein.clear_hydration();
        auto grid = protein.get_grid();

        auto vol = grid::exv::RawGridExv::create(grid);
        auto vol2 = grid::exv::RawGridWithSurfaceExv::create(grid, false);

        REQUIRE(vol.interior.size() == vol2.interior.size());
        REQUIRE(vol.surface.size() == vol2.surface.size());

        for (unsigned int i = 0; i < vol.interior.size(); i++) {
            CHECK(vol.interior[i] == vol2.interior[i]);
        }
    }
}