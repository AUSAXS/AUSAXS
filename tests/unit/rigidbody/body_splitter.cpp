#include <catch2/catch_test_macros.hpp>

#include <rigidbody/BodySplitter.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <io/File.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::data;

TEST_CASE("BodySplitter::split by chain") {
    SECTION("split multi-chain protein") {
        auto molecule = BodySplitter::split(io::File("tests/files/SASDJG5.pdb"));
        CHECK(molecule.get_bodies().size() == 2);
    }

    SECTION("split single-chain protein") {
        auto molecule = BodySplitter::split(io::File("tests/files/2epe.pdb"));
        CHECK(molecule.get_bodies().size() == 1);
    }
}

TEST_CASE("BodySplitter::split by indices") {
    SECTION("split at single index") {
        auto molecule = BodySplitter::split(io::File("tests/files/2epe.pdb"), {50});
        CHECK(molecule.get_bodies().size() == 2);
    }

    SECTION("split at multiple indices") {
        auto molecule = BodySplitter::split(io::File("tests/files/2epe.pdb"), {30, 60, 90});
        CHECK(molecule.get_bodies().size() == 4);
    }

    SECTION("split creates non-empty bodies") {
        auto molecule = BodySplitter::split(io::File("tests/files/2epe.pdb"), {50});
        CHECK(molecule.get_body(0).size_atom() > 0);
        CHECK(molecule.get_body(1).size_atom() > 0);
    }
}
