#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactorInfo.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("form_factor::to_string") {
    SECTION("to_string and from_string round trip") {
        for (int i = 0; i < total_ff_count; ++i) {
            auto type = static_cast<form_factor_t>(i);
            CHECK(from_string(to_string(type)) == type);
        }
    }

    SECTION("untabulated types") {
        CHECK(to_string(form_factor_t::COUNT) == "CNT");
        CHECK(to_string(form_factor_t::UNKNOWN) == "UNK");
    }
}

TEST_CASE("form_factor::get_info") {
    SECTION("the vacuum form factor of the element at q=0 matches its number of electrons") {
        for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
            auto type = static_cast<form_factor_t>(i);
            const auto& info = get_info(type);
            const auto& c = info.xray_coefficients;
            double f0 = c.c;
            for (double a : c.a) {f0 += a;}
            INFO("form factor " << to_string(type));
            CHECK_THAT(f0, Catch::Matchers::WithinAbs(constants::charge::nuclear::get_charge(info.element), 0.01));
        }
    }

    SECTION("invalid types throw") {
        CHECK_THROWS(get_info(form_factor_t::COUNT));
        CHECK_THROWS(get_info(form_factor_t::UNKNOWN));
    }
}

TEST_CASE("form_factor::get_type") {
    SECTION("each bare element maps back to its own type") {
        for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
            auto type = static_cast<form_factor_t>(i);
            const auto& info = get_info(type);
            if (type == form_factor_t::OTHER || info.geometry.hydrogens != 0) {continue;}
            CHECK(get_type(info.element) == type);
        }
    }

    SECTION("elements without their own type map to OTHER") {
        CHECK(get_type(constants::atom_t::Ar) == form_factor_t::OTHER);
        CHECK(get_type(constants::atom_t::Au) == form_factor_t::OTHER);
        CHECK(get_type(constants::atom_t::unknown) == form_factor_t::OTHER);
    }

    SECTION("atomic group takes priority") {
        CHECK(get_type(constants::atom_t::C, constants::atomic_group_t::CH) == form_factor_t::CH);
        CHECK(get_type(constants::atom_t::C, constants::atomic_group_t::CH2) == form_factor_t::CH2);
        CHECK(get_type(constants::atom_t::C, constants::atomic_group_t::CH3) == form_factor_t::CH3);
        CHECK(get_type(constants::atom_t::N, constants::atomic_group_t::NH) == form_factor_t::NH);
        CHECK(get_type(constants::atom_t::N, constants::atomic_group_t::NH2) == form_factor_t::NH2);
        CHECK(get_type(constants::atom_t::N, constants::atomic_group_t::NH3) == form_factor_t::NH3);
        CHECK(get_type(constants::atom_t::O, constants::atomic_group_t::OH) == form_factor_t::OH);
        CHECK(get_type(constants::atom_t::S, constants::atomic_group_t::SH) == form_factor_t::SH);
    }

    SECTION("unknown atomic group falls back to the element") {
        CHECK(get_type(constants::atom_t::H, constants::atomic_group_t::unknown) == form_factor_t::H);
        CHECK(get_type(constants::atom_t::C, constants::atomic_group_t::unknown) == form_factor_t::C);
    }
}

TEST_CASE("form_factor::to_atom_type") {
    SECTION("types without an element map to argon") {
        CHECK(to_atom_type(form_factor_t::EXCLUDED_VOLUME) == constants::atom_t::Ar);
        CHECK(to_atom_type(form_factor_t::COUNT) == constants::atom_t::Ar);
        CHECK(to_atom_type(form_factor_t::UNKNOWN) == constants::atom_t::Ar);
    }
}

TEST_CASE("constants::mass::get_mass") {
    CHECK(constants::mass::get_mass(form_factor_t::COUNT) == 0);
}

TEST_CASE("constants::radius::get_vdw_radius") {
    SECTION("unknown types occupy the space of a carbon") {
        CHECK(constants::radius::get_vdw_radius(form_factor_t::UNKNOWN) == constants::radius::get_vdw_radius(constants::atom_t::C));
    }

    SECTION("types without an element throw") {
        CHECK_THROWS(constants::radius::get_vdw_radius(form_factor_t::EXCLUDED_VOLUME));
        CHECK_THROWS(constants::radius::get_vdw_radius(form_factor_t::COUNT));
    }
}
