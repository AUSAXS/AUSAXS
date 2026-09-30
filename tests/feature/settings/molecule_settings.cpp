#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <settings/All.h>

#include <support/temp_file.h>

using namespace ausaxs;

TEST_CASE("MoleculeSettings::allow_unknown_residues") {
	SECTION("true") {
        settings::molecule::allow_unknown_residues = true;
        REQUIRE_NOTHROW(data::Molecule("tests/files/diamond.pdb"));
	}

	SECTION("false") {
        settings::molecule::allow_unknown_residues = false;
        REQUIRE_THROWS(data::Molecule("tests/files/diamond.pdb"));
	}
}
