#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <data/Molecule.h>
#include <io/detail/structure/XYZReader.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace data;

TEST_CASE("XYZReader::read") {
    SECTION("simple") {
        // Create a temporary XYZ file
        io::File tmp("temp/tests/io/temp.xyz");
        tmp.create(
            "3\n"
            "Comment line\n"
            "H 0.0 0.0 0.0\n"
            "O 0.0 0.0 1.0\n"
            "H 1.0 0.0 0.0\n"
            "Au 1.2 3.4 5.6\n"
        );

        auto structure = io::detail::xyz::read(tmp);
        REQUIRE(structure.atoms.size() == 4);

        REQUIRE(structure.atoms[0].element == constants::atom_t::H);
        REQUIRE(structure.atoms[0].coordinates().x() == 0.0);
        REQUIRE(structure.atoms[0].coordinates().y() == 0.0);
        REQUIRE(structure.atoms[0].coordinates().z() == 0.0);

        REQUIRE(structure.atoms[1].element == constants::atom_t::O);
        REQUIRE(structure.atoms[1].coordinates().x() == 0.0);
        REQUIRE(structure.atoms[1].coordinates().y() == 0.0);
        REQUIRE(structure.atoms[1].coordinates().z() == 1.0);

        REQUIRE(structure.atoms[2].element == constants::atom_t::H);
        REQUIRE(structure.atoms[2].coordinates().x() == 1.0);
        REQUIRE(structure.atoms[2].coordinates().y() == 0.0);
        REQUIRE(structure.atoms[2].coordinates().z() == 0.0);

        REQUIRE(structure.atoms[3].element == constants::atom_t::Au);
        REQUIRE(structure.atoms[3].coordinates().x() == 1.2);
        REQUIRE(structure.atoms[3].coordinates().y() == 3.4);
        REQUIRE(structure.atoms[3].coordinates().z() == 5.6);
    }

    SECTION("real file") {
        io::File file("tests/files/Au.xyz");
        auto structure = io::detail::xyz::read(file);
        REQUIRE(structure.atoms.size() == 201);

        REQUIRE(structure.atoms[0].element == constants::atom_t::Au);
        REQUIRE_THAT(structure.atoms[0].coordinates().x(), Catch::Matchers::WithinAbs(-3.81837662, 1e-6));
        REQUIRE_THAT(structure.atoms[0].coordinates().y(), Catch::Matchers::WithinAbs(-7.63675324, 1e-6));
        REQUIRE_THAT(structure.atoms[0].coordinates().z(), Catch::Matchers::WithinAbs(0., 1e-6));

        REQUIRE(structure.atoms[1].element == constants::atom_t::Au);
        REQUIRE_THAT(structure.atoms[1].coordinates().x(), Catch::Matchers::WithinAbs(3.81837662, 1e-6));
        REQUIRE_THAT(structure.atoms[1].coordinates().y(), Catch::Matchers::WithinAbs(-7.63675324, 1e-6));
        REQUIRE_THAT(structure.atoms[1].coordinates().z(), Catch::Matchers::WithinAbs(0., 1e-6));
    }
}