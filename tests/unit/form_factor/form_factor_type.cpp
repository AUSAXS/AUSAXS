#include <catch2/catch_test_macros.hpp>

#include <form_factor/FormFactorType.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("form_factor_t::total_ff_count") {
    SECTION("count value") {
        int count = total_ff_count;
        CHECK(count == static_cast<int>(form_factor_t::COUNT));
        CHECK(count > 0);
    }
}

TEST_CASE("form_factor_t::bins") {
    SECTION("exv_bin") {
        CHECK(exv_bin == static_cast<int>(form_factor_t::EXCLUDED_VOLUME));
    }

    SECTION("water_bin") {
        CHECK(water_bin == static_cast<int>(form_factor_t::OH));
    }
}
