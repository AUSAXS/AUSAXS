// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/FormFactorManager.h>

#include <constants/ConstantsAxes.h>
#include <data/Body.h>  // IWYU pragma: keep
#include <data/Molecule.h>
#include <form_factor/FormFactorConcepts.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/detail/LookupHelpers.h>
#include <utility/Logging.h>

#include <algorithm>
#include <cassert>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    std::unique_ptr<manager::detail::ActiveTables> active_tables;
    using ff_profile_t = std::array<double, constants::axes::q_axis.bins>; // A single form factor evaluated over the default q axis.
    using profile_set_t = std::array<ff_profile_t, form_factor::total_ff_count>; // One such profile per active form factor slot.

    /**
     * @brief Evaluate every active atomic form factor over the default q axis.
     *        Evaluating first & then multiplying the results is faster than evaluating each product individually. 
     */
    template<FormFactorLookupType FormFactorLookup>
    profile_set_t evaluate_atomic_profiles(const std::array<int, form_factor::total_ff_count>& ff_indices) {
        profile_set_t profiles{};
        for (int i = 0; i < form_factor::get_active_count(); ++i) {
            const auto& ff = FormFactorLookup::get(static_cast<form_factor_t>(ff_indices[i]));
            for (int q = 0; q < static_cast<int>(constants::axes::q_axis.bins); ++q) {
                profiles[i][q] = ff.evaluate(constants::axes::q_vals[q]);
            }
        }
        return profiles;
    }

    /**
     * @brief Evaluate every explicit excluded volume form factor of the current EXV set over the default q axis.
     */
    profile_set_t evaluate_exv_profiles(const std::array<int, form_factor::total_ff_count>& ff_indices) {
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        profile_set_t profiles{};
        for (int i = start_index_for_explicit_exv(); i < form_factor::get_active_count(); ++i) {
            auto ff = exv_set.get(static_cast<form_factor_t>(ff_indices[i]));
            for (int q = 0; q < static_cast<int>(constants::axes::q_axis.bins); ++q) {
                profiles[i][q] = ff.evaluate(constants::axes::q_vals[q]);
            }
        }
        return profiles;
    }

    /**
     * @brief Generate an atomic form factor product table.
     */
    lookup::table_t generate_atomic_table(const profile_set_t& atomic) {
        lookup::table_t table;
        for (int i = 0; i < form_factor::get_active_count(); ++i) {
            for (int j = 0; j < i; ++j) {
                table.index(i, j) = FormFactorProduct(atomic[i], atomic[j]);
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(atomic[i], atomic[i]);
        }
        return table;
    }

    /**
     * @brief Generate an excluded volume form factor product table (exv-exv). This is a symmetric table.
     */
    lookup::table_t generate_exv_table(const profile_set_t& exv) {
        lookup::table_t table;
        for (int i = start_index_for_explicit_exv(); i < form_factor::get_active_count(); ++i) {
            for (int j = start_index_for_explicit_exv(); j < i; ++j) {
                table.index(i, j) = FormFactorProduct(exv[i], exv[j]);
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(exv[i], exv[i]);
        }
        return table;
    }

    /**
     * @brief Generate a cross form factor product table (atomic-exv).
     */
    lookup::table_t generate_cross_table(const profile_set_t& atomic, const profile_set_t& exv) {
        lookup::table_t table;
        for (int i = 0; i < form_factor::get_active_count(); ++i) {
            for (int j = start_index_for_explicit_exv(); j < form_factor::get_active_count(); ++j) {
                table.index(i, j) = FormFactorProduct(atomic[i], exv[j]);
            }
        }
        return table;
    }
}

manager::detail::ActiveTables::ActiveTables(const std::array<int, form_factor::total_ff_count>& ff_indices, int active_count) 
    : active_count(active_count), ff_indices(ff_indices)
{
    // must come first; the profile evaluations and table generators below only fill the active sub-block, which they read from here
    form_factor::detail::active_ff_count = active_count;

    const auto raw_profiles        = evaluate_atomic_profiles<lookup::detail::RawFormFactorLookup>(this->ff_indices);
    const auto normalized_profiles = evaluate_atomic_profiles<lookup::detail::NormalizedFormFactorLookup>(this->ff_indices);
    const auto exv_profiles        = evaluate_exv_profiles(this->ff_indices);

    this->raw_atomic_table         = generate_atomic_table(raw_profiles);
    this->raw_cross_table          = generate_cross_table(raw_profiles, exv_profiles);
    this->raw_exv_table            = generate_exv_table(exv_profiles);
    this->normalized_atomic_table  = generate_atomic_table(normalized_profiles);
    this->normalized_cross_table   = generate_cross_table(normalized_profiles, exv_profiles);
}

observer_ptr<const manager::detail::ActiveTables> manager::get_active_product_tables() noexcept {
    if (!active_tables) { // initialize default tables
        std::array<int, form_factor::total_ff_count> default_indices;
        std::iota(default_indices.begin(), default_indices.end(), 0);
        active_tables = std::make_unique<detail::ActiveTables>(default_indices, form_factor::total_ff_count);
    }
    return active_tables.get();
}

std::vector<int> manager::get_active_mapping() {
    auto ff_indices = get_active_product_tables()->ff_indices;
    std::vector<int> mapping(form_factor::total_ff_count, -1);
    for (int i = 0; i < form_factor::get_active_count(); ++i) {
        mapping[ff_indices[i]] = i;
    }

    // form factors not in the active set fall back to the OTHER slot
    int other_slot = mapping[static_cast<int>(form_factor::form_factor_t::OTHER)];
    assert(other_slot != -1 && "OTHER must always be part of the active form factor set.");
    for (auto& m : mapping) {if (m == -1) {m = other_slot;}}
    return mapping;
}

void manager::detail::use_form_factors(std::vector<int> ff_indices) {
    assert(!ff_indices.empty() && "Custom form factors cannot be empty.");
    assert(ff_indices.size() <= form_factor::total_ff_count && "Custom form factors cannot exceed the total number of available form factors.");

    // ensure form_factor_t::OTHER is always present
    constexpr int other = static_cast<int>(form_factor::form_factor_t::OTHER);
    if (std::ranges::find(ff_indices, other) == ff_indices.end()) {
        assert(ff_indices.size() < form_factor::total_ff_count && "Cannot append OTHER to a full form factor set.");
        ff_indices.push_back(other);
    }

    std::array<int, form_factor::total_ff_count> ff_indices_array;
    std::ranges::copy(ff_indices, ff_indices_array.begin());
    std::fill(ff_indices_array.begin() + ff_indices.size(), ff_indices_array.end(), static_cast<int>(form_factor::form_factor_t::OTHER));
    active_tables = std::make_unique<detail::ActiveTables>(ff_indices_array, ff_indices.size());
}

void manager::use_form_factors(data::Molecule& molecule) {
    std::vector<int> ff_counts(form_factor::total_ff_count, 0);
    for (auto& a : molecule.iterate_atoms()) {
        ++ff_counts[static_cast<int>(a.form_factor_type())];
    }
    // ensure excluded volume and water are always at the front of the list, and OTHER is always at the end, regardless of abundance
    ff_counts[static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME)] = std::numeric_limits<int>::max();
    ff_counts[static_cast<int>(form_factor::form_factor_t::WATER)] = std::numeric_limits<int>::max()-1;
    ff_counts[static_cast<int>(form_factor::form_factor_t::OTHER)] = std::numeric_limits<int>::min();

    std::vector<int> ff_indices(form_factor::total_ff_count);
    std::iota(ff_indices.begin(), ff_indices.end(), 0);
    std::ranges::sort(ff_indices, [&ff_counts](int a, int b) {return ff_counts[a] > ff_counts[b];});

    // Truncate to the form factors actually present. The sort above places EXCLUDED_VOLUME and WATER first (forced), then every type with a non-zero atom 
    // count in descending order, then the absent types, and finally OTHER. Everything from the first absent type onwards is dead weight and therefore removed. 
    int n_present = 0;
    for (int i = 2; i < static_cast<int>(ff_indices.size()); ++i) {
        if (ff_counts[ff_indices[i]] <= 0) {break;}
        ++n_present;
    }
    ff_indices.resize(std::min<int>(2 + n_present + 1, form_factor::total_ff_count));
    ff_indices.back() = static_cast<int>(form_factor::form_factor_t::OTHER); // OTHER will never be selected, so it is safe to assign it here

    if (logging::logging_enabled()) {
        std::string log_msg = "Setting form factors based on detected molecular composition:";
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::EXCLUDED_VOLUME) + " (forced)";
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::WATER) + " (forced)";
        for (int i = 2; i < static_cast<int>(ff_indices.size())-1; ++i) {
            log_msg += "\n\t" + form_factor::to_string(static_cast<form_factor_t>(ff_indices[i])) + " with count " + std::to_string(ff_counts[ff_indices[i]]);
        }
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::OTHER) + " (forced)";
        logging::log(log_msg);
    }

    detail::use_form_factors(std::move(ff_indices));
}

void manager::rebuild() {
    if (!active_tables) {return;} // lazy init will pick up the new EXV set
    active_tables = std::make_unique<detail::ActiveTables>(
        std::array<int, form_factor::total_ff_count>(active_tables->ff_indices),
        active_tables->active_count
    );
}