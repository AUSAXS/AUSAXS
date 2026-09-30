#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <data/Molecule.h>
#include <io/detail/structure/XYZReader.h>
#include <settings/All.h>

#include <support/temp_file.h>

using namespace ausaxs;
using namespace data;

TEST_CASE("XYZReader::read") {
    SECTION("simple") {
        test::TempFile tmp(".xyz",
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

        auto reduced = structure.reduced_representation();
        for (int i = 0; i < static_cast<int>(reduced.atoms.size()); ++i) {
            REQUIRE_THAT(
                static_cast<double>(reduced.atoms[i].weight()),
                Catch::Matchers::WithinRel(structure.atoms[i].effective_charge, 1e-6)
            );
        }
        // gold must not be weighted as an unknown form factor would be
        REQUIRE(reduced.atoms[3].weight() > 70);
    }

    SECTION("real file") {
        io::File file("tests/files/carbon_sphere.xyz");
        auto structure = io::detail::xyz::read(file);
        REQUIRE(structure.atoms.size() == 1985);

        REQUIRE(structure.atoms[0].element == constants::atom_t::C);
        REQUIRE_THAT(structure.atoms[0].coordinates().x(), Catch::Matchers::WithinAbs(-18.35212500, 1e-6));
        REQUIRE_THAT(structure.atoms[0].coordinates().y(), Catch::Matchers::WithinAbs(-6.11737500, 1e-6));
        REQUIRE_THAT(structure.atoms[0].coordinates().z(), Catch::Matchers::WithinAbs(-4.07825000, 1e-6));

        REQUIRE(structure.atoms[1].element == constants::atom_t::C);
        REQUIRE_THAT(structure.atoms[1].coordinates().x(), Catch::Matchers::WithinAbs(-18.35212500, 1e-6));
        REQUIRE_THAT(structure.atoms[1].coordinates().y(), Catch::Matchers::WithinAbs(-4.07825000, 1e-6));
        REQUIRE_THAT(structure.atoms[1].coordinates().z(), Catch::Matchers::WithinAbs(-6.11737500, 1e-6));
    }
}

TEST_CASE("XYZReader::read does not disable implicit hydrogens for later structures", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = true;

    Molecule before("tests/files/2epe.pdb");
    double expected = before.get_total_atomic_charge();

    {   // the .xyz itself must still get no implicit hydrogens
        io::File file("tests/files/carbon_sphere.xyz");
        auto structure = io::detail::xyz::read(file);
        REQUIRE_FALSE(structure.supports_implicit_hydrogens);
        structure.add_implicit_hydrogens();
        for (const auto& a : structure.atoms) {REQUIRE(a.element == constants::atom_t::C);}
    }

    Molecule sphere("tests/files/carbon_sphere.xyz");
    REQUIRE(settings::molecule::implicit_hydrogens);

    Molecule after("tests/files/2epe.pdb");
    REQUIRE_THAT(after.get_total_atomic_charge(), Catch::Matchers::WithinRel(expected, 1e-12));
}