#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/atoms/Atom.h>
#include <settings/All.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>

using namespace ausaxs;

TEST_CASE("test_unknown_form_factor: UNKNOWN form factors with Fraser excluded volume model") {
    // Set settings BEFORE creating the molecule so the histogram manager uses Fraser method
    settings::fit::fit_hydration = true;
    settings::fit::fit_excluded_volume = true;
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser;

    // Create atoms without form factor information (like molecule_from_arrays does)
    // Using data::Atom creates atoms without form factor information, which will be converted
    // to UNKNOWN form factors when added to a Body
    std::vector<data::Atom> atoms;
    atoms.emplace_back(Vector3<double>{0.0, 0.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{1.0, 0.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{0.0, 1.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{0.0, 0.0, 1.0}, 1.0);
    atoms.emplace_back(Vector3<double>{-1.0, 0.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{0.0, -1.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{0.0, 0.0, -1.0}, 1.0);
    atoms.emplace_back(Vector3<double>{1.0, 1.0, 1.0}, 1.0);
    atoms.emplace_back(Vector3<double>{-1.0, -1.0, -1.0}, 1.0);

    REQUIRE(atoms.size() == 9);

    // Create molecule from atoms (this will convert them to AtomFF with UNKNOWN form factors)
    auto molecule = data::Molecule({data::Body{atoms}});
    REQUIRE(molecule.size_atom() == 9);

    // Verify that atoms have UNKNOWN form factors (when no form factor info is available)
    auto& body = molecule.get_bodies()[0];
    auto& atoms_ff = body.get_atoms();
    CHECK(atoms_ff[0].form_factor_type() == form_factor::form_factor_t::UNKNOWN);

    // This should throw an informative error when trying to use UNKNOWN form factors with Fraser
    REQUIRE_THROWS_WITH([&]() {
        auto hist = molecule.get_histogram();
        auto I = hist->debye_transform();
    }(), Catch::Matchers::ContainsSubstring("UNKNOWN form factor"));
}

TEST_CASE("test_unknown_form_factor_simple: UNKNOWN form factors with Simple excluded volume model") {
    // Set Simple ExV method - this should work fine with UNKNOWN form factors
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;

    // Create atoms without form factor information
    std::vector<data::Atom> atoms;
    atoms.emplace_back(Vector3<double>{0.0, 0.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{1.0, 0.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{0.0, 1.0, 0.0}, 1.0);
    atoms.emplace_back(Vector3<double>{0.0, 0.0, 1.0}, 1.0);

    // Create molecule
    auto molecule = data::Molecule({data::Body{atoms}});
    REQUIRE(molecule.size_atom() == 4);

    // Verify atoms have UNKNOWN form factors
    auto& body = molecule.get_bodies()[0];
    auto& atoms_ff = body.get_atoms();
    CHECK(atoms_ff[0].form_factor_type() == form_factor::form_factor_t::UNKNOWN);

    // This should work fine with Simple ExV model (doesn't need form factor info)
    REQUIRE_NOTHROW([&]() {
        auto hist = molecule.get_histogram();
        auto I = hist->debye_transform();
        REQUIRE(I.size() > 0);
    }());
}
