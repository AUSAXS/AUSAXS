#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include <utility/Settings.h>

TEST_CASE("io") {
    SECTION("write_settings") {
        setting::write("temp/settings/settings.txt");
    }

    SECTION("read_settings") {
        setting::read("temp/settings/settings.txt");
    }
}