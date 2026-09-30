#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <settings/All.h>

#include <support/temp_file.h>

using namespace ausaxs;

TEST_CASE("settings::write and settings::read") {
    SECTION("write_settings") {
        test::TempFile path(".txt");
        settings::write(path);
    }

    SECTION("read_settings") {
        test::TempFile path(".txt");
        settings::write(path);
        settings::read(path);
    }
}
