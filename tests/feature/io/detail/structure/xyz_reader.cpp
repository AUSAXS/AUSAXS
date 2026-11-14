#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <data/Body.h>
#include <data/Molecule.h>
#include <io/detail/structure/XYZReader.h>
#include <settings/All.h>

#include <support/temp_file.h>

using namespace ausaxs;
using namespace data;

TEST_CASE("XYZReader::read") {
    SECTION("simple") {
        test::TempFile tmp(".xyz",
            "4\n"
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

        data::Body body(structure);
        for (int i = 0; i < body.size_atom(); ++i) {
            REQUIRE_THAT(
                static_cast<double>(body.get_atom(i).weight()),
                Catch::Matchers::WithinRel(constants::charge::get_ff_charge(body.get_atom(i).form_factor_type(), structure.atoms[i].element), 1e-6)
            );
        }
        // gold must not be weighted as an unknown form factor would be
        REQUIRE(body.get_atom(3).weight() > 70);
    }

    SECTION("comment line is never an atom") {
        auto comment = GENERATE(
            std::string("Au 1 2 3"),
            std::string("Ag Au nanoparticle model"),
            std::string("{'Ag': 1, 'Au': 1}"),
            std::string("Properties=species:S:1:pos:R:3 pbc=\"F F F\"")
        );
        test::TempFile tmp(".xyz", "2\n" + comment + "\nAg 0 0 0\nAu 2.9 0 0\n");
        auto structure = io::detail::xyz::read(tmp);
        REQUIRE(structure.atoms.size() == 2);
        CHECK(structure.atoms[0].element == constants::atom_t::Ag);
        CHECK(structure.atoms[1].element == constants::atom_t::Au);
        CHECK(structure.atoms[1].coordinates().x() == 2.9);
    }

    SECTION("extended xyz with extra per-atom columns") {
        test::TempFile tmp(".xyz",
            "2\n"
            "Lattice=\"20 0 0 0 20 0 0 0 20\" Properties=species:S:1:pos:R:3:tags:I:1 pbc=\"T T T\"\n"
            "Au 10.0 10.0 10.0 0\n"
            "Au 12.0 10.0 10.0 1\n"
        );
        auto structure = io::detail::xyz::read(tmp);
        REQUIRE(structure.atoms.size() == 2);
        CHECK(structure.atoms[1].coordinates().x() == 12.0);
    }

    SECTION("malformed files are refused") {
        auto content = GENERATE(
            std::string(""),                                // no atom count
            std::string("0\ncomment\n"),                    // no atoms
            std::string("C 0 0 0\nC 1 0 0\n"),               // no header
            std::string("3\ncomment\nC 0 0 0\nC 1 0 0\n"),  // fewer atoms than the count
            std::string("1\ncomment\nC 0 0 0\nC 1 0 0\n"),  // more atoms than the count
            std::string("1\nframe 1\nC 0 0 0\n1\nframe 2\nC 1 0 0\n"), // multi-frame
            std::string("1\ncomment\nC 0 0\n"),             // missing coordinate
            std::string("1\ncomment\nXx 0 0 0\n")           // unknown element
        );
        test::TempFile tmp(".xyz", content);
        REQUIRE_THROWS_AS(io::detail::xyz::read(tmp), except::io_error);
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

TEST_CASE("XYZReader::read does not disable implicit hydrogens for later structures") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = true;

    Molecule before("tests/files/2epe.pdb");
    double expected = before.get_total_atomic_charge();

    {   // the .xyz itself must still get no implicit hydrogens
        io::File file("tests/files/carbon_sphere.xyz");
        auto structure = io::detail::xyz::read(file);
        REQUIRE_FALSE(structure.supports_implicit_hydrogens);
        for (const auto& a : structure.reduced_representation().atoms) {REQUIRE(a.form_factor_type() == form_factor::form_factor_t::C);}
    }

    Molecule sphere("tests/files/carbon_sphere.xyz");
    REQUIRE(settings::molecule::implicit_hydrogens);

    Molecule after("tests/files/2epe.pdb");
    REQUIRE_THAT(after.get_total_atomic_charge(), Catch::Matchers::WithinRel(expected, 1e-12));
}