#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/SymmetryManager.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/Symmetry.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("SymmetryManager::calculate") {
    settings::molecule::implicit_hydrogens = false;
    settings::general::threads = 1;

    SECTION("simple") {
        record::Atom a(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);
        a.set_effective_charge(1);
        Molecule m({a});
        m.get_body(0).add_symmetry({Vector3<double>(1, 0, 0), Vector3<double>(0, 0, 0)});

        hist::detail::SymmetryManager sm;
        auto h = sm.calculate<false>(m)->get_total_counts();

        REQUIRE(h[0] == 2);
        REQUIRE(h[std::round(1*constants::axes::d_inv_width)] == 2);
    }
}