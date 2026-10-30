#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <constants/Constants.h>
#include <constants/ConstantsProperties.h>
#include <constants/ConstantsSI.h>
#include <constants/ConstantsUnits.h>

using namespace ausaxs;

TEST_CASE("units: the two conversion namespaces agree") {
    // 1 Da expressed in kg, reached through either namespace
    CHECK_THAT(constants::unit::external::gm*constants::SI::mass::gm, Catch::Matchers::WithinRel(constants::SI::mass::u, 1e-15));
    CHECK_THAT(constants::unit::external::mg*constants::SI::mass::mg, Catch::Matchers::WithinRel(constants::SI::mass::u, 1e-15));

    // 1 A expressed in m, reached through either namespace
    CHECK_THAT(
        constants::unit::external::cm*constants::SI::length::cm,
        Catch::Matchers::WithinRel(constants::SI::length::A, 1e-15)
    );
    CHECK_THAT(
        constants::unit::external::nm*constants::SI::length::nm,
        Catch::Matchers::WithinRel(constants::SI::length::A, 1e-15)
    );
    CHECK_THAT(
        constants::unit::external::mL*constants::SI::volume::cm3,
        Catch::Matchers::WithinRel(constants::SI::volume::A3, 1e-15)
    );
}

TEST_CASE("units: densities are in internal units") {
    auto to_gcm3 = [] (double rho) {return rho*constants::unit::external::gm/constants::unit::external::mL;};

    // water is 0.9982071 g/cm^3 and 0.334 e/A^3; a Da is ~1.66e-24 g and an A^3 is 1e-24 cm^3,
    // so a density in Da/A^3 is a factor ~1.66 below the same density in g/cm^3.
    CHECK_THAT(constants::mass::density::water,   Catch::Matchers::WithinRel(0.6011344, 1e-6));
    CHECK_THAT(constants::mass::density::protein, Catch::Matchers::WithinRel(0.8129890, 1e-6));
    CHECK_THAT(to_gcm3(constants::mass::density::water),   Catch::Matchers::WithinRel(0.9982071, 1e-12));
    CHECK_THAT(to_gcm3(constants::mass::density::protein), Catch::Matchers::WithinRel(1.35,      1e-12));

    // a protein is denser than water, and both are of order 1 in internal units - an SI density
    // would be three orders of magnitude out
    CHECK(constants::mass::density::water < constants::mass::density::protein);
    CHECK(constants::mass::density::protein < 10);
    CHECK(0.1 < constants::mass::density::water);
}
