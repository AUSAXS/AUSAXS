#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <constants/ConstantsAxes.h>
#include <data/Body.h>  // IWYU pragma: keep
#include <data/Molecule.h>
#include <form_factor/FormFactor.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <settings/All.h>
#include <support/form_factor_helper.h>

#include <algorithm>
#include <concepts>
#include <limits>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::form_factor;

static const std::vector<int>& identity() {
    static std::vector<int> identity;
    if (identity.empty()) {
        identity = std::vector<int>(total_ff_count);
        std::iota(identity.begin(), identity.end(), 0);
    }
    return identity;
}

TEST_CASE("form_factor_manager: full identity selection") {
    settings::exv::exv_method = settings::exv::ExvMethod::Simple; // the Fraser-based models remove the types without an excluded volume
    auto original_max = settings::form_factor::max_types;
    settings::form_factor::max_types = total_ff_count;
    manager::detail::use_form_factors(identity());
    const auto* tables = manager::get_active_product_tables();
    REQUIRE(tables != nullptr);

    SECTION("active_count equals total_ff_count") {
        REQUIRE(tables->active_count == form_factor::total_ff_count);
        REQUIRE(get_active_count() == form_factor::total_ff_count);
    }

    SECTION("ff_indices are identity") {
        for (int i = 0; i < form_factor::total_ff_count; ++i) {
            REQUIRE(tables->ff_indices[i] == static_cast<int>(i));
        }
    }

    SECTION("identity mapping") {
        auto mapping = manager::get_active_mapping();
        REQUIRE(mapping.size() == total_ff_count);
        for (int i = 0; i < total_ff_count; ++i) {
            REQUIRE(mapping[i] == static_cast<int>(i));
        }
    }

    settings::form_factor::max_types = original_max;
}

TEST_CASE("form_factor::get_active_count") {
    test::form_factor::use_random_form_factors();
    REQUIRE(get_active_count() == manager::get_active_product_tables()->active_count);

    SECTION("reflects custom subset") {
        manager::detail::use_form_factors({
            static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
            static_cast<int>(form_factor_t::WATER),
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::OTHER)
        });
        REQUIRE(get_active_count() == 4);
        REQUIRE(get_active_count() == manager::get_active_product_tables()->active_count);
    }

    SECTION("OTHER is appended to a subset that omits it") {
        // every inactive type is folded onto the OTHER slot, so it must always be part of the active set
        manager::detail::use_form_factors({
            static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
            static_cast<int>(form_factor_t::WATER),
            static_cast<int>(form_factor_t::C)
        });
        REQUIRE(get_active_count() == 4);
        REQUIRE(manager::get_active_product_tables()->ff_indices[3] == static_cast<int>(form_factor_t::OTHER));
    }
}

TEST_CASE("form_factor_manager::get_active_mapping custom subset") {
    const int exv   = static_cast<int>(form_factor_t::EXCLUDED_VOLUME);
    const int water = static_cast<int>(form_factor_t::WATER);
    const int C     = static_cast<int>(form_factor_t::C);
    const int other = static_cast<int>(form_factor_t::OTHER);

    manager::detail::use_form_factors({exv, water, C});

    auto mapping = manager::get_active_mapping();
    REQUIRE(mapping.size() == total_ff_count);

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

    SECTION("OTHER slot is the last active slot") {
        // OTHER is appended to the requested set, so it takes the last slot of the active prefix.
        // The trailing padding is deliberately *not* mapped: those slots lie outside the histograms,
        // which are sized to the active count.
        REQUIRE(get_active_count() == 4);
        REQUIRE(mapping[other] == 3);
    }
}

TEST_CASE("form_factor_manager::detail::use_form_factors padding") {
    manager::detail::use_form_factors({
        static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
        static_cast<int>(form_factor_t::WATER),
        static_cast<int>(form_factor_t::C),
        static_cast<int>(form_factor_t::N)
    });

    const auto* tables = manager::get_active_product_tables();

    SECTION("active_count reflects explicit count plus the appended OTHER") {
        REQUIRE(tables->active_count == 5);
    }

    SECTION("explicit slots are preserved") {
        REQUIRE(tables->ff_indices[0] == static_cast<int>(form_factor_t::EXCLUDED_VOLUME));
        REQUIRE(tables->ff_indices[1] == static_cast<int>(form_factor_t::WATER));
        REQUIRE(tables->ff_indices[2] == static_cast<int>(form_factor_t::C));
        REQUIRE(tables->ff_indices[3] == static_cast<int>(form_factor_t::N));
    }

    SECTION("OTHER is appended after the explicit slots") {
        REQUIRE(tables->ff_indices[4] == static_cast<int>(form_factor_t::OTHER));
    }

    SECTION("trailing slots are padded with OTHER") {
        for (int i = 5; i < form_factor::total_ff_count; ++i) {
            REQUIRE(tables->ff_indices[i] == static_cast<int>(form_factor_t::OTHER));
        }
    }
}

TEST_CASE("form_factor_manager::use_form_factors(Molecule) ordering") {
    data::Molecule molecule("tests/files/2epe.pdb");
    manager::use_form_factors(molecule);

    const auto* tables = manager::get_active_product_tables();

    SECTION("EXV is always slot 0") {
        REQUIRE(tables->ff_indices[0] == static_cast<int>(form_factor_t::EXCLUDED_VOLUME));
    }

    SECTION("WATER is always slot 1") {
        REQUIRE(tables->ff_indices[1] == static_cast<int>(form_factor_t::WATER));
    }

    SECTION("OTHER is always the last active slot") {
        REQUIRE(tables->ff_indices[tables->active_count - 1] == static_cast<int>(form_factor_t::OTHER));
    }

    SECTION("active set is truncated to the types actually present") {
        // the active set is [EXV, WATER, <present types>, OTHER]; every absent type is dropped
        std::vector<int> counts(static_cast<int>(form_factor_t::UNKNOWN) + 1, 0);
        for (const auto& a : molecule.iterate_atoms()) {
            ++counts[static_cast<int>(a.form_factor_type())];
        }

        int expected = 3; // EXV and WATER are forced to the front, OTHER to the back
        for (int t = 0; t < static_cast<int>(total_ff_count); ++t) {
            if (t == static_cast<int>(form_factor_t::EXCLUDED_VOLUME)) {continue;}
            if (t == static_cast<int>(form_factor_t::WATER)) {continue;}
            if (t == static_cast<int>(form_factor_t::OTHER)) {continue;}
            if (counts[t] > 0) {++expected;}
        }

        REQUIRE(tables->active_count == expected);
    }

    SECTION("every type present in the molecule has its own slot") {
        auto mapping = manager::get_active_mapping();
        int other_slot = mapping[static_cast<int>(form_factor_t::OTHER)];
        for (const auto& a : molecule.iterate_atoms()) {
            int t = static_cast<int>(a.form_factor_type());
            if (t == static_cast<int>(form_factor_t::OTHER)) {continue;} // OTHER is always folded onto its own slot
            REQUIRE(mapping[t] != other_slot);
        }
    }
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
    int count_before = manager::get_active_product_tables()->active_count;

    // capture one exv table value before rebuild; the exv table is only filled from start_index_for_explicit_exv() onwards
    double exv_val_before = manager::get_active_product_tables()->raw_exv_table.index(start_index_for_explicit_exv(), start_index_for_explicit_exv()).evaluate(0);

    manager::rebuild();

    const auto* tables = manager::get_active_product_tables();

    SECTION("ff_indices unchanged after rebuild") {
        REQUIRE(tables->ff_indices == indices_before);
    }

    SECTION("active_count unchanged after rebuild") {
        REQUIRE(tables->active_count == count_before);
    }

    SECTION("table values reproduced identically after rebuild with same EXV set") {
        double exv_val_after = tables->raw_exv_table.index(start_index_for_explicit_exv(), start_index_for_explicit_exv()).evaluate(0);
        REQUIRE(exv_val_before == exv_val_after);
    }
}

TEST_CASE("form_factor_manager::rebuild after EXV set change updates exv table") {
    test::form_factor::use_random_form_factors();

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

TEST_CASE("form_factor_manager: product tables hold the product their indices name") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models

    auto check = [] <std::invocable<form_factor_t, double> Row, std::invocable<form_factor_t, double> Col> (
        const lookup::table_t& table, int i0, int j0, const Row& row, const Col& col
    ) {
        const auto* tables = manager::get_active_product_tables();
        for (int i = i0; i < get_active_count(); ++i) {
            for (int j = j0; j < get_active_count(); ++j) {
                auto ti = static_cast<form_factor_t>(tables->ff_indices[i]);
                auto tj = static_cast<form_factor_t>(tables->ff_indices[j]);
                for (int q = 0; q < static_cast<int>(constants::axes::q_axis.bins); ++q) {
                    double expected = row(ti, constants::axes::q_vals[q])*col(tj, constants::axes::q_vals[q]);
                    REQUIRE_THAT(table.index(i, j).evaluate(q), Catch::Matchers::WithinRel(expected, 1e-12));
                }
            }
        }
    };

    auto raw = [] (form_factor_t t, double q) {return lookup::atomic::raw::get(t).evaluate(q);};
    auto normalized = [] (form_factor_t t, double q) {return lookup::atomic::normalized::get(t).evaluate(q);};
    auto exv = [] (form_factor_t t, double q) {return ExvTableManager::get_current_exv_form_factor_set().get(t).evaluate(q);};
    int s0 = start_index_for_explicit_exv();

    SECTION("random set") {
        test::form_factor::use_random_form_factors();
        const auto* tables = manager::get_active_product_tables();
        check(tables->raw_atomic_table,        0,  0,  raw,        raw);
        check(tables->normalized_atomic_table, 0,  0,  normalized, normalized);
        check(tables->raw_cross_table,         0,  s0, raw,        exv);
        check(tables->normalized_cross_table,  0,  s0, normalized, exv);
        check(tables->raw_exv_table,           s0, s0, exv,        exv);
    }

    SECTION("truncated set from a molecule") {
        data::Molecule molecule("tests/files/2epe.pdb");
        manager::use_form_factors(molecule);
        REQUIRE(get_active_count() < total_ff_count); // the truncation has to actually bite
        const auto* tables = manager::get_active_product_tables();
        check(tables->raw_atomic_table,        0,  0,  raw,        raw);
        check(tables->normalized_atomic_table, 0,  0,  normalized, normalized);
        check(tables->raw_cross_table,         0,  s0, raw,        exv);
        check(tables->normalized_cross_table,  0,  s0, normalized, exv);
        check(tables->raw_exv_table,           s0, s0, exv,        exv);
    }
}

TEST_CASE("form_factor_manager: Fraser only uses form factors with an excluded volume") {
    const int CH3   = static_cast<int>(form_factor_t::CH3);
    const int other = static_cast<int>(form_factor_t::OTHER);
    auto original_method = settings::exv::exv_method.value;
    auto original_set = settings::exv::exv_set.value;
    auto original_max = settings::form_factor::max_types;
    auto original_fraction = settings::form_factor::min_fraction;
    settings::form_factor::max_types = total_ff_count;
    manager::detail::use_form_factors(identity());

    // a volume set without CH3
    auto set = constants::exv::MinimumFluctuation_implicit_H;
    set.volumes[CH3].reset();
    ExvTableManager::set_custom_exv_table(set);

    // not every type has a volume in the base set either, so count those that remain
    int with_volume = 1; // EXCLUDED_VOLUME is always kept
    for (int i = start_index_for_explicit_exv(); i < total_ff_count; ++i) {
        if (set.contains(static_cast<form_factor_t>(i))) {++with_volume;}
    }

    auto is_active = [] (int type) {
        const auto* tables = manager::get_active_product_tables();
        return std::ranges::find(tables->ff_indices.begin(), tables->ff_indices.begin() + tables->active_count, type) != tables->ff_indices.begin() + tables->active_count;
    };

    SECTION("Fraser removes the type, and maps it onto OTHER") {
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
        CHECK(get_active_count() == with_volume);
        CHECK_FALSE(is_active(CH3));
        CHECK(is_active(other));
        auto mapping = manager::get_active_mapping();
        CHECK(mapping[CH3] == mapping[other]);
    }

    SECTION("other models keep the type") {
        settings::exv::exv_method = settings::exv::ExvMethod::Grid;
        CHECK(get_active_count() == total_ff_count);
        CHECK(is_active(CH3));
    }

    SECTION("switching the model restores the requested selection") {
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
        CHECK_FALSE(is_active(CH3));
        settings::exv::exv_method = settings::exv::ExvMethod::Grid;
        CHECK(is_active(CH3));
    }

    SECTION("molecule-based selection skips the type") {
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
        data::Molecule molecule("tests/files/2epe.pdb");
        auto atoms = molecule.iterate_atoms();
        REQUIRE(std::ranges::any_of(atoms, [] (const auto& a) {return a.form_factor_type() == form_factor_t::CH3;}));
        manager::use_form_factors(molecule);
        CHECK_FALSE(is_active(CH3));
    }

    SECTION("the type does not take up a slot of the molecule-based selection") {
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
        data::Molecule molecule("tests/files/2epe.pdb");
        std::vector<int> counts(total_ff_count, 0);
        for (const auto& a : molecule.iterate_atoms()) {
            if (form_factor::detail::is_tabulated(a.form_factor_type())) {++counts[static_cast<int>(a.form_factor_type())];}
        }

        // limit the slots to exactly the types at least as abundant as CH3, so CH3 would take the last one if it were not skipped
        int at_least_ch3 = 0, below_ch3 = 0;
        for (int t = start_index_for_explicit_exv(); t < total_ff_count; ++t) {
            if (t == static_cast<int>(form_factor_t::WATER) || t == other || counts[t] == 0) {continue;}
            if (counts[CH3] <= counts[t]) {++at_least_ch3;}
            else {++below_ch3;}
        }
        REQUIRE(0 < below_ch3); // some type must be left to take the freed slot
        settings::form_factor::max_types = 3 + at_least_ch3;
        settings::form_factor::min_fraction = 0;

        manager::use_form_factors(molecule);
        CHECK_FALSE(is_active(CH3));
        CHECK(get_active_count() == settings::form_factor::max_types);
    }

    settings::exv::exv_method = original_method;
    settings::exv::exv_set = original_set;
    settings::form_factor::max_types = original_max;
    settings::form_factor::min_fraction = original_fraction;
}

TEST_CASE("form_factor_manager: max_types is a hard limit") {
    const int exv   = static_cast<int>(form_factor_t::EXCLUDED_VOLUME);
    const int water = static_cast<int>(form_factor_t::WATER);
    const int C     = static_cast<int>(form_factor_t::C);
    const int N     = static_cast<int>(form_factor_t::N);
    const int O     = static_cast<int>(form_factor_t::O);
    const int other = static_cast<int>(form_factor_t::OTHER);
    auto original_max = settings::form_factor::max_types;
    settings::form_factor::max_types = 5;

    SECTION("a selection within the limit is accepted") {
        manager::detail::use_form_factors({exv, water, C, N, other});
        CHECK(get_active_count() == 5);
    }

    SECTION("a selection above the limit is rejected") {
        CHECK_THROWS(manager::detail::use_form_factors({exv, water, C, N, O, other}));
    }

    SECTION("the appended OTHER counts towards the limit") {
        CHECK_THROWS(manager::detail::use_form_factors({exv, water, C, N, O}));
    }

    settings::form_factor::max_types = original_max;
}

TEST_CASE("form_factor_manager: use_form_factors(Molecule) folds excess and rare types onto OTHER") {
    settings::exv::exv_method = settings::exv::ExvMethod::Simple; // the Fraser-based models remove the types without an excluded volume
    auto original_max = settings::form_factor::max_types;
    auto original_fraction = settings::form_factor::min_fraction;
    data::Molecule molecule("tests/files/2epe.pdb");

    // the atom counts of the present types, excluding the forced EXCLUDED_VOLUME, WATER, and OTHER
    std::vector<int> counts(total_ff_count, 0);
    for (const auto& a : molecule.iterate_atoms()) {
        if (form_factor::detail::is_tabulated(a.form_factor_type())) {++counts[static_cast<int>(a.form_factor_type())];}
    }
    std::vector<int> present;
    for (int t = 0; t < total_ff_count; ++t) {
        if (t == static_cast<int>(form_factor_t::EXCLUDED_VOLUME) || t == static_cast<int>(form_factor_t::WATER) || t == static_cast<int>(form_factor_t::OTHER)) {continue;}
        if (0 < counts[t]) {present.push_back(t);}
    }
    REQUIRE(4 < present.size());

    auto has_own_slot = [] (int type) {
        auto mapping = manager::get_active_mapping();
        return mapping[type] != mapping[static_cast<int>(form_factor_t::OTHER)];
    };

    SECTION("the slot limit keeps the most abundant types") {
        settings::form_factor::max_types = 6;
        settings::form_factor::min_fraction = 0;
        manager::use_form_factors(molecule);
        const auto* tables = manager::get_active_product_tables();
        REQUIRE(tables->active_count == 6);
        CHECK(tables->ff_indices[0] == static_cast<int>(form_factor_t::EXCLUDED_VOLUME));
        CHECK(tables->ff_indices[1] == static_cast<int>(form_factor_t::WATER));
        CHECK(tables->ff_indices[5] == static_cast<int>(form_factor_t::OTHER));

        int min_kept = std::numeric_limits<int>::max(), max_folded = 0;
        for (int t : present) {
            if (has_own_slot(t)) {min_kept = std::min(min_kept, counts[t]);}
            else {max_folded = std::max(max_folded, counts[t]);}
        }
        CHECK(max_folded <= min_kept);
    }

    SECTION("types below the minimum fraction are folded onto OTHER") {
        settings::form_factor::max_types = total_ff_count;
        int rarest = *std::ranges::min_element(present, {}, [&counts] (int t) {return counts[t];});
        settings::form_factor::min_fraction = (counts[rarest] + 0.5)/static_cast<double>(molecule.size_atom());
        manager::use_form_factors(molecule);
        for (int t : present) {
            CHECK(has_own_slot(t) == (counts[rarest] < counts[t]));
        }
    }

    settings::form_factor::max_types = original_max;
    settings::form_factor::min_fraction = original_fraction;
}
