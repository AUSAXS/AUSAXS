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

    SECTION("one copy") {
        record::Atom a(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);
        a.set_effective_charge(1);
        Molecule m({a});
        m.get_body(0).add_symmetry({Vector3<double>(1, 0, 0), Vector3<double>(0, 0, 0)});

        hist::detail::SymmetryManager sm;
        auto h = sm.calculate<false>(m)->get_total_counts();

        int bin = std::round(1*constants::axes::d_inv_width);
        REQUIRE(bin < static_cast<int>(h.size()));
        REQUIRE(h[0] == 2);
        for (int i = 1; i < bin; ++i) {
            REQUIRE(h[i] == 0);
        }
        REQUIRE(h[std::round(1*constants::axes::d_inv_width)] == 2);
        for (int i = bin+1; i < static_cast<int>(h.size()); ++i) {
            REQUIRE(h[i] == 0);
        }
    }

    SECTION("two copies") {
        record::Atom a(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);
        a.set_effective_charge(1);
        Molecule m({a});
        m.get_body(0).add_symmetry({Vector3<double>(-1, 0, 0), Vector3<double>(0, 0, 0)});
        m.get_body(0).add_symmetry({Vector3<double>( 1, 0, 0), Vector3<double>(0, 0, 0)});

        hist::detail::SymmetryManager sm;
        auto h = sm.calculate<false>(m)->get_total_counts();

        int bin1 = std::round(1*constants::axes::d_inv_width);
        int bin2 = std::round(2*constants::axes::d_inv_width);
        REQUIRE(bin2 < static_cast<int>(h.size()));
        REQUIRE(h[0] == 3);
        for (int i = 1; i < bin1; ++i) {
            REQUIRE(h[i] == 0);
        }
        REQUIRE(h[bin1] == 4);
        for (int i = bin1+1; i < bin2; ++i) {
            REQUIRE(h[i] == 0);
        }
        REQUIRE(h[bin2] == 2);
        for (int i = bin2+1; i < static_cast<int>(h.size()); ++i) {
            REQUIRE(h[i] == 0);
        }
    }

    SECTION("four copies") {
        record::Atom a(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);
        a.set_effective_charge(1);
        Molecule m({a});
        m.get_body(0).add_symmetry({Vector3<double>(-1, 0, 0), Vector3<double>(0, 0, 0)});
        m.get_body(0).add_symmetry({Vector3<double>( 1, 0, 0), Vector3<double>(0, 0, 0)});
        m.get_body(0).add_symmetry({Vector3<double>( 0,-1, 0), Vector3<double>(0, 0, 0)});
        m.get_body(0).add_symmetry({Vector3<double>( 0, 1, 0), Vector3<double>(0, 0, 0)});

        hist::detail::SymmetryManager sm;
        auto h = sm.calculate<false>(m)->get_total_counts();

        int bin1 = std::round(1*constants::axes::d_inv_width);
        int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
        int bin3 = std::round(2*constants::axes::d_inv_width);
        REQUIRE(bin3 < static_cast<int>(h.size()));
        REQUIRE(h[0] == 5);
        for (int i = 1; i < bin1; ++i) {
            REQUIRE(h[i] == 0);
        }
        REQUIRE(h[bin1] == 8);
        for (int i = bin1+1; i < bin2; ++i) {
            REQUIRE(h[i] == 0);
        }
        REQUIRE(h[bin2] == 8);
        for (int i = bin2+1; i < bin3; ++i) {
            REQUIRE(h[i] == 0);
        }
        REQUIRE(h[bin3] == 4);
        for (int i = bin3+1; i < static_cast<int>(h.size()); ++i) {
            REQUIRE(h[i] == 0);
        }
    }
}