#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <data/Atom.h>
#include <data/Water.h>
#include <settings/All.h>
#include <utility/Constants.h>

TEST_CASE("Water::create_new_water") {
    SECTION("water") {
        Water w1 = Water::create_new_water(Vector3<double>({1, 2, 3}));
        CHECK(w1.serial == -1);
        CHECK(w1.name == "O");
        CHECK(w1.altLoc == "");
        CHECK(w1.resName == "HOH");
        CHECK(w1.chainID == "");
        CHECK(w1.resSeq == -1);
        CHECK(w1.iCode == "");
        CHECK(w1.coords == Vector3({1, 2, 3}));
        CHECK(w1.occupancy == 1);
        CHECK(w1.tempFactor == 0);
        CHECK(w1.element == "O");
        CHECK(w1.charge == "");
        CHECK(w1.get_type() == RecordType::WATER);
        CHECK(w1.is_water() == true);
    }
}