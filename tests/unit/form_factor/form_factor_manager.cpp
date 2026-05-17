#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/FormFactorType.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

#include <numeric>

using namespace ausaxs;
using namespace ausaxs::form_factor;

const std::vector<int>& identity() {
    static std::vector<int> identity;
    if (identity.empty()) {
        identity = std::vector<int>(get_total_ff_count());
        std::iota(identity.begin(), identity.end(), 0);
    }
    return identity;
}

TEST_CASE("form_factor_manager::get_active_product_tables lazy init") {
    auto* tables = manager::get_active_product_tables();
    REQUIRE(tables != nullptr);

    SECTION("active_count equals max_ff_types for identity set") {
        REQUIRE(tables->active_count == static_cast<unsigned int>(settings::form_factor::max_ff_types));
    }

    SECTION("ff_indices are identity") {
        for (unsigned int i = 0; i < static_cast<unsigned int>(settings::form_factor::max_ff_types); ++i) {
            REQUIRE(tables->ff_indices[i] == static_cast<int>(i));
        }
    }
}

TEST_CASE("form_factor::get_active_count") {
    REQUIRE(get_active_count() == manager::get_active_product_tables()->active_count);

    SECTION("reflects custom subset") {
        manager::detail::use_form_factors({
            static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
            static_cast<int>(form_factor_t::WATER),
            static_cast<int>(form_factor_t::C)
        });
        REQUIRE(get_active_count() == 3);
        REQUIRE(get_active_count() == manager::get_active_product_tables()->active_count);
        manager::detail::use_form_factors(identity());
    }
}

TEST_CASE("form_factor_manager::get_active_mapping default") {
    auto mapping = manager::get_active_mapping();
    REQUIRE(mapping.size() == get_total_ff_count());

    SECTION("identity mapping") {
        for (unsigned int i = 0; i < get_total_ff_count(); ++i) {
            REQUIRE(mapping[i] == static_cast<int>(i));
        }
    }
}

TEST_CASE("form_factor_manager::get_active_mapping custom subset") {
    const int exv   = static_cast<int>(form_factor_t::EXCLUDED_VOLUME);
    const int water = static_cast<int>(form_factor_t::WATER);
    const int C     = static_cast<int>(form_factor_t::C);
    const int other = static_cast<int>(form_factor_t::OTHER);

    manager::detail::use_form_factors({exv, water, C});

    auto mapping = manager::get_active_mapping();
    REQUIRE(mapping.size() == get_total_ff_count());

    SECTION("active types map to correct slots") {
        REQUIRE(mapping[exv]   == 0);
        REQUIRE(mapping[water] == 1);
        REQUIRE(mapping[C]     == 2);
    }

    SECTION("inactive types fall back to the OTHER slot") {
        // types not in the active set must map to a real, in-bounds slot (the OTHER slot)
        REQUIRE(mapping[static_cast<int>(form_factor_t::N)] == mapping[other]);
        REQUIRE(mapping[static_cast<int>(form_factor_t::H)] == mapping[other]);
    }

    SECTION("OTHER slot is last in the padded array") {
        // padding fills trailing slots with OTHER, so the mapping for OTHER is overwritten repeatedly and ends at the last padded index
        REQUIRE(mapping[other] == settings::form_factor::max_ff_types - 1);
    }

    manager::detail::use_form_factors(identity());
}

TEST_CASE("form_factor_manager::detail::use_form_factors padding") {
    manager::detail::use_form_factors({
        static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
        static_cast<int>(form_factor_t::WATER),
        static_cast<int>(form_factor_t::C),
        static_cast<int>(form_factor_t::N)
    });

    auto* tables = manager::get_active_product_tables();

    SECTION("active_count reflects explicit count") {
        REQUIRE(tables->active_count == 4);
    }

    SECTION("explicit slots are preserved") {
        REQUIRE(tables->ff_indices[0] == static_cast<int>(form_factor_t::EXCLUDED_VOLUME));
        REQUIRE(tables->ff_indices[1] == static_cast<int>(form_factor_t::WATER));
        REQUIRE(tables->ff_indices[2] == static_cast<int>(form_factor_t::C));
        REQUIRE(tables->ff_indices[3] == static_cast<int>(form_factor_t::N));
    }

    SECTION("trailing slots are padded with OTHER") {
        for (int i = 4; i < settings::form_factor::max_ff_types; ++i) {
            REQUIRE(tables->ff_indices[i] == static_cast<int>(form_factor_t::OTHER));
        }
    }

    manager::detail::use_form_factors(identity());
}

TEST_CASE("form_factor_manager::use_form_factors(Molecule) ordering") {
    data::Molecule molecule("tests/files/2epe.pdb");
    manager::use_form_factors(molecule);

    auto* tables = manager::get_active_product_tables();

    SECTION("EXV is always slot 0") {
        REQUIRE(tables->ff_indices[0] == static_cast<int>(form_factor_t::EXCLUDED_VOLUME));
    }

    SECTION("WATER is always slot 1") {
        REQUIRE(tables->ff_indices[1] == static_cast<int>(form_factor_t::WATER));
    }

    SECTION("OTHER is always the last slot") {
        REQUIRE(tables->ff_indices[settings::form_factor::max_ff_types - 1] == static_cast<int>(form_factor_t::OTHER));
    }

    SECTION("active_count equals max_ff_types for a real molecule with diverse atoms") {
        // 2epe.pdb contains many different atom types, expect all slots populated
        REQUIRE(tables->active_count == static_cast<unsigned int>(settings::form_factor::max_ff_types));
    }

    manager::detail::use_form_factors(identity());
}

TEST_CASE("form_factor_manager::rebuild preserves indices and regenerates tables") {
    // set a custom subset
    manager::detail::use_form_factors({
        static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
        static_cast<int>(form_factor_t::WATER),
        static_cast<int>(form_factor_t::C),
        static_cast<int>(form_factor_t::N)
    });

    auto indices_before = manager::get_active_product_tables()->ff_indices;
    unsigned int count_before = manager::get_active_product_tables()->active_count;

    // capture one exv table value before rebuild
    double exv_val_before = manager::get_active_product_tables()->raw_exv_table.index(0, 0).evaluate(0);

    manager::rebuild();

    auto* tables = manager::get_active_product_tables();

    SECTION("ff_indices unchanged after rebuild") {
        REQUIRE(tables->ff_indices == indices_before);
    }

    SECTION("active_count unchanged after rebuild") {
        REQUIRE(tables->active_count == count_before);
    }

    SECTION("table values reproduced identically after rebuild with same EXV set") {
        double exv_val_after = tables->raw_exv_table.index(0, 0).evaluate(0);
        REQUIRE(exv_val_before == exv_val_after);
    }

    manager::detail::use_form_factors(identity());
}

TEST_CASE("form_factor_manager::rebuild after EXV set change updates exv table") {
    // capture exv product at (1,1) which is WATER vs WATER excluded volume — should differ between sets
    double exv_val_default = manager::get_active_product_tables()->raw_exv_table.index(1, 1).evaluate(10);

    // directly update the setting value to avoid the stale-read bug in the callback, then manually call rebuild so it reads the new value
    settings::exv::exv_set.value = settings::exv::ExvSet::Traube;
    manager::rebuild();

    double exv_val_traube = manager::get_active_product_tables()->raw_exv_table.index(1, 1).evaluate(10);

    REQUIRE(exv_val_default != exv_val_traube);

    settings::exv::exv_set.value = settings::exv::ExvSet::Default;
    manager::rebuild();
}