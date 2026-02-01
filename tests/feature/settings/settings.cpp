#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <utility/Utility.h>
#include <data/Molecule.h>
#include <constants/ConstantsAxes.h>
#include <settings/All.h>

using namespace ausaxs;

TEST_CASE("settings") {
    SECTION("write_settings") {
        settings::write("temp/settings/settings.txt");
    }

    SECTION("read_settings") {
        settings::read("temp/settings/settings.txt");
    }
}

TEST_CASE("allow_unknown_residues") {
	SECTION("true") {
        settings::molecule::allow_unknown_residues = true;
        REQUIRE_NOTHROW(data::Molecule("tests/files/diamond.pdb"));
	}

	SECTION("false") {
        settings::molecule::allow_unknown_residues = false;
        REQUIRE_THROWS(data::Molecule("tests/files/diamond.pdb"));
	}
}