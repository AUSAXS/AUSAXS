// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/detail/FormFactorProductBase.h>
#include <form_factor/lookup/detail/LookupHelpers.h>
#include <settings/FormFactorSettings.h>
#include <container/Container2D.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/Logging.h>

#include <cassert>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    std::unique_ptr<manager::detail::ActiveTables> active_tables;
}

manager::detail::ActiveTables::ActiveTables(std::array<int, settings::form_factor::max_ff_types>&& ff_indices, unsigned int active_count) 
    : ff_indices(std::move(ff_indices)), active_count(active_count)
{
    this->raw_atomic_table         = lookup::detail::generate_atomic_table<lookup::detail::RawFormFactorLookup>(this->ff_indices);
    this->raw_cross_table          = lookup::detail::generate_cross_table<lookup::detail::RawFormFactorLookup>(this->ff_indices);
    this->raw_exv_table            = lookup::detail::generate_exv_table(this->ff_indices);
    this->normalized_atomic_table  = lookup::detail::generate_atomic_table<lookup::detail::NormalizedFormFactorLookup>(this->ff_indices);
    this->normalized_cross_table   = lookup::detail::generate_cross_table<lookup::detail::NormalizedFormFactorLookup>(this->ff_indices);
}

observer_ptr<const manager::detail::ActiveTables> manager::get_active_product_tables() noexcept {
    if (!active_tables) { // initialize default tables
        std::array<int, settings::form_factor::max_ff_types> default_indices;
        std::iota(default_indices.begin(), default_indices.end(), 0);
        active_tables = std::make_unique<detail::ActiveTables>(std::move(default_indices), settings::form_factor::max_ff_types);
    }
    return active_tables.get();
}

std::vector<int> manager::get_active_mapping() {
    auto ff_indices = get_active_product_tables()->ff_indices;
    std::vector<int> mapping(form_factor::get_total_ff_count(), settings::form_factor::max_ff_types);
    assert(mapping.size() == ff_indices.size() && "Mapping size should match total form factor count.");
    for (unsigned int i = 0; i < ff_indices.size(); ++i) {
        mapping[ff_indices[i]] = i;
    }
    return mapping;
}

unsigned int form_factor::get_active_count() noexcept {
    return manager::get_active_product_tables()->active_count;
}

void manager::detail::use_form_factors(std::vector<int> ff_indices) {
    assert(!ff_indices.empty() && "Custom form factors cannot be empty.");
    assert(ff_indices.size() <= form_factor::get_total_ff_count() && "Custom form factors cannot exceed the total number of available form factors.");

    std::array<int, settings::form_factor::max_ff_types> ff_indices_array;
    std::copy(ff_indices.begin(), ff_indices.end(), ff_indices_array.begin());
    std::fill(ff_indices_array.begin() + ff_indices.size(), ff_indices_array.end(), static_cast<int>(form_factor::form_factor_t::OTHER));
    active_tables = std::make_unique<detail::ActiveTables>(std::move(ff_indices_array), ff_indices.size());
}

void manager::use_form_factors(data::Molecule& molecule) {
    std::vector<int> ff_counts(form_factor::get_total_ff_count(), 0);
    for (auto& a : molecule.iterate_atoms()) {
        ++ff_counts[static_cast<int>(a.form_factor_type())];
    }
    // ensure excluded volume and water are always at the front of the list, and OTHER is always at the end, regardless of abundance
    ff_counts[static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME)] = std::numeric_limits<int>::max();
    ff_counts[static_cast<int>(form_factor::form_factor_t::WATER)] = std::numeric_limits<int>::max()-1;
    ff_counts[static_cast<int>(form_factor::form_factor_t::OTHER)] = std::numeric_limits<int>::min();

    std::vector<int> ff_indices(form_factor::get_total_ff_count());
    std::iota(ff_indices.begin(), ff_indices.end(), 0);
    std::sort(ff_indices.begin(), ff_indices.end(), [&ff_counts](int a, int b) {return ff_counts[a] > ff_counts[b];});
    ff_indices.resize(settings::form_factor::max_ff_types);
    ff_indices.back() = static_cast<int>(form_factor::form_factor_t::OTHER); // OTHER will never be selected, so it is safe to assign it here

    if (logging::logging_enabled()) {
        std::string log_msg = "Setting form factors based on detected molecular composition:";
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::EXCLUDED_VOLUME) + " (forced)";
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::WATER) + " (forced)";
        for (unsigned int i = 2; i < ff_indices.size()-1; ++i) {
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
        std::array<int, settings::form_factor::max_ff_types>(active_tables->ff_indices),
        active_tables->active_count
    );
}