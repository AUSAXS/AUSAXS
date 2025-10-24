#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactorType.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("form_factor_t::to_string") {
    SECTION("basic elements") {
        CHECK(to_string(form_factor_t::H) == "H");
        CHECK(to_string(form_factor_t::C) == "C");
        CHECK(to_string(form_factor_t::N) == "N");
        CHECK(to_string(form_factor_t::O) == "O");
        CHECK(to_string(form_factor_t::S) == "S");
    }

    SECTION("atomic groups") {
        CHECK(to_string(form_factor_t::CH) == "CH");
        CHECK(to_string(form_factor_t::CH2) == "CH2");
        CHECK(to_string(form_factor_t::CH3) == "CH3");
        CHECK(to_string(form_factor_t::NH) == "NH");
        CHECK(to_string(form_factor_t::NH2) == "NH2");
        CHECK(to_string(form_factor_t::NH3) == "NH3");
        CHECK(to_string(form_factor_t::OH) == "OH");
        CHECK(to_string(form_factor_t::SH) == "SH");
    }

    SECTION("special types") {
        CHECK(to_string(form_factor_t::OTHER) == "OTH");
        CHECK(to_string(form_factor_t::EXCLUDED_VOLUME) == "EXV");
        CHECK(to_string(form_factor_t::COUNT) == "CNT");
        CHECK(to_string(form_factor_t::UNKNOWN) == "UNK");
    }
}

TEST_CASE("form_factor_t::get_count") {
    SECTION("count value") {
        unsigned int count = get_count();
        CHECK(count == static_cast<unsigned int>(form_factor_t::COUNT));
        CHECK(count > 0);
    }
}

TEST_CASE("form_factor_t::get_count_without_excluded_volume") {
    SECTION("count without excluded volume") {
        unsigned int count = get_count_without_excluded_volume();
        CHECK(count == get_count() - 1);
        CHECK(count == static_cast<unsigned int>(form_factor_t::COUNT) - 1);
    }
}

TEST_CASE("form_factor_t::get_type") {
    SECTION("from atom_t") {
        CHECK(get_type(constants::atom_t::H) == form_factor_t::H);
        CHECK(get_type(constants::atom_t::C) == form_factor_t::C);
        CHECK(get_type(constants::atom_t::N) == form_factor_t::N);
        CHECK(get_type(constants::atom_t::O) == form_factor_t::O);
        CHECK(get_type(constants::atom_t::S) == form_factor_t::S);
        CHECK(get_type(constants::atom_t::Ar) == form_factor_t::OTHER);
    }

    SECTION("from atom_t and atomic_group_t") {
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

        SECTION("fallback to atom_t") {
            CHECK(get_type(constants::atom_t::H, constants::atomic_group_t::unknown) == form_factor_t::H);
            CHECK(get_type(constants::atom_t::C, constants::atomic_group_t::unknown) == form_factor_t::C);
        }
    }
}

TEST_CASE("form_factor_t::to_atom_type") {
    SECTION("basic elements") {
        CHECK(to_atom_type(form_factor_t::H) == constants::atom_t::H);
        CHECK(to_atom_type(form_factor_t::C) == constants::atom_t::C);
        CHECK(to_atom_type(form_factor_t::N) == constants::atom_t::N);
        CHECK(to_atom_type(form_factor_t::O) == constants::atom_t::O);
        CHECK(to_atom_type(form_factor_t::S) == constants::atom_t::S);
    }

    SECTION("carbon groups") {
        CHECK(to_atom_type(form_factor_t::CH) == constants::atom_t::C);
        CHECK(to_atom_type(form_factor_t::CH2) == constants::atom_t::C);
        CHECK(to_atom_type(form_factor_t::CH3) == constants::atom_t::C);
    }

    SECTION("nitrogen groups") {
        CHECK(to_atom_type(form_factor_t::NH) == constants::atom_t::N);
        CHECK(to_atom_type(form_factor_t::NH2) == constants::atom_t::N);
        CHECK(to_atom_type(form_factor_t::NH3) == constants::atom_t::N);
    }

    SECTION("oxygen groups") {
        CHECK(to_atom_type(form_factor_t::OH) == constants::atom_t::O);
    }

    SECTION("sulfur groups") {
        CHECK(to_atom_type(form_factor_t::SH) == constants::atom_t::S);
    }

    SECTION("other") {
        CHECK(to_atom_type(form_factor_t::OTHER) == constants::atom_t::Ar);
    }
}

TEST_CASE("constants::mass::get_mass") {
    SECTION("basic elements") {
        CHECK(constants::mass::get_mass(form_factor_t::H) > 0);
        CHECK(constants::mass::get_mass(form_factor_t::C) > 0);
        CHECK(constants::mass::get_mass(form_factor_t::N) > 0);
        CHECK(constants::mass::get_mass(form_factor_t::O) > 0);
        CHECK(constants::mass::get_mass(form_factor_t::S) > 0);
    }

    SECTION("atomic groups have greater mass") {
        CHECK(constants::mass::get_mass(form_factor_t::CH) > constants::mass::get_mass(form_factor_t::C));
        CHECK(constants::mass::get_mass(form_factor_t::CH2) > constants::mass::get_mass(form_factor_t::CH));
        CHECK(constants::mass::get_mass(form_factor_t::CH3) > constants::mass::get_mass(form_factor_t::CH2));
        
        CHECK(constants::mass::get_mass(form_factor_t::NH) > constants::mass::get_mass(form_factor_t::N));
        CHECK(constants::mass::get_mass(form_factor_t::NH2) > constants::mass::get_mass(form_factor_t::NH));
        CHECK(constants::mass::get_mass(form_factor_t::NH3) > constants::mass::get_mass(form_factor_t::NH2));
        
        CHECK(constants::mass::get_mass(form_factor_t::OH) > constants::mass::get_mass(form_factor_t::O));
        CHECK(constants::mass::get_mass(form_factor_t::SH) > constants::mass::get_mass(form_factor_t::S));
    }

    SECTION("special types") {
        CHECK(constants::mass::get_mass(form_factor_t::OTHER) > 0);
        CHECK(constants::mass::get_mass(form_factor_t::EXCLUDED_VOLUME) == 0);
        CHECK(constants::mass::get_mass(form_factor_t::COUNT) == 0);
    }
}

TEST_CASE("constants::radius::get_vdw_radius") {
    SECTION("basic elements") {
        CHECK(constants::radius::get_vdw_radius(form_factor_t::H) > 0);
        CHECK(constants::radius::get_vdw_radius(form_factor_t::C) > 0);
        CHECK(constants::radius::get_vdw_radius(form_factor_t::N) > 0);
        CHECK(constants::radius::get_vdw_radius(form_factor_t::O) > 0);
        CHECK(constants::radius::get_vdw_radius(form_factor_t::S) > 0);
    }

    SECTION("atomic groups inherit base atom radius") {
        CHECK(constants::radius::get_vdw_radius(form_factor_t::CH) == constants::radius::get_vdw_radius(form_factor_t::C));
        CHECK(constants::radius::get_vdw_radius(form_factor_t::CH2) == constants::radius::get_vdw_radius(form_factor_t::C));
        CHECK(constants::radius::get_vdw_radius(form_factor_t::CH3) == constants::radius::get_vdw_radius(form_factor_t::C));
        
        CHECK(constants::radius::get_vdw_radius(form_factor_t::NH) == constants::radius::get_vdw_radius(form_factor_t::N));
        CHECK(constants::radius::get_vdw_radius(form_factor_t::NH2) == constants::radius::get_vdw_radius(form_factor_t::N));
        CHECK(constants::radius::get_vdw_radius(form_factor_t::NH3) == constants::radius::get_vdw_radius(form_factor_t::N));
        
        CHECK(constants::radius::get_vdw_radius(form_factor_t::OH) == constants::radius::get_vdw_radius(form_factor_t::O));
        CHECK(constants::radius::get_vdw_radius(form_factor_t::SH) == constants::radius::get_vdw_radius(form_factor_t::S));
    }

    SECTION("special types") {
        CHECK(constants::radius::get_vdw_radius(form_factor_t::OTHER) > 0);
        CHECK(constants::radius::get_vdw_radius(form_factor_t::UNKNOWN) == 0);
    }
}

TEST_CASE("constants::charge::nuclear::get_charge") {
    SECTION("basic elements") {
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::H) == 1);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::C) == 6);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::N) == 7);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::O) == 8);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::S) == 16);
    }

    SECTION("atomic groups have combined charge") {
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::CH) == 7);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::CH2) == 8);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::CH3) == 9);
        
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::NH) == 8);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::NH2) == 9);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::NH3) == 10);
        
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::OH) == 9);
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::SH) == 17);
    }

    SECTION("other") {
        CHECK(constants::charge::nuclear::get_charge(form_factor_t::OTHER) == 18);
    }
}

TEST_CASE("form_factor_t::bins") {
    SECTION("exv_bin") {
        CHECK(exv_bin == static_cast<int>(form_factor_t::EXCLUDED_VOLUME));
    }

    SECTION("water_bin") {
        CHECK(water_bin == static_cast<int>(form_factor_t::OH));
    }
}
