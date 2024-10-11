#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <settings/All.h>

TEST_CASE("settings") {
    SECTION("write_settings") {
        settings::write("temp/settings/settings.txt");
    }

    SECTION("read_settings") {
        settings::read("temp/settings/settings.txt");
    }
}