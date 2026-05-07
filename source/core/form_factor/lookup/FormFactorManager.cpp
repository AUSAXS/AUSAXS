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

observer_ptr<const FormFactorManager::_CustomTables> FormFactorManager::get_custom_tables() noexcept {
    return custom_tables.get();
}

#if CONSTEXPR_TABLES
    #define CONST constexpr
#else
    #define CONST const
#endif
namespace {
    // Default tables. These _must_ be initialized statically to ensure default settings. 
    CONST auto raw_atomic_table        = lookup::detail::generate_atomic_table<lookup::detail::RawFormFactorLookup>();
    CONST auto raw_exv_table           = lookup::detail::generate_exv_table(true);
    CONST auto raw_cross_table         = lookup::detail::generate_cross_table<lookup::detail::RawFormFactorLookup>(true);
    CONST auto normalized_atomic_table = lookup::detail::generate_atomic_table<lookup::detail::NormalizedFormFactorLookup>();
    CONST auto normalized_cross_table  = lookup::detail::generate_cross_table<lookup::detail::NormalizedFormFactorLookup>(true);
}

std::vector<int> FormFactorManager::get_ff_mapping() {
    // default to OTHER which is always last
    auto ff_indices = is_using_custom_form_factors() ? custom_tables->ff_indices : get_ff_indices();
    unsigned int count = is_using_custom_form_factors() ? custom_tables->active_count : ff_indices.size();
    std::vector<int> mapping(form_factor::get_total_ff_count(), settings::form_factor::max_ff_types);
    for (unsigned int i = 0; i < count; ++i) {
        mapping[ff_indices[i]] = i;
    }
    return mapping;
}

unsigned int FormFactorManager::get_active_count() noexcept {
    if (is_using_custom_form_factors() && custom_tables) {
        return custom_tables->active_count;
    }
    return settings::form_factor::max_ff_types;
}

const lookup::atomic::table_t& FormFactorManager::normalized_atomic_table() noexcept {
    if (is_using_custom_form_factors()) {
        refresh_custom_state();
        return custom_tables->custom_normalized_atomic_table;
    }
    return ::normalized_atomic_table;
}

const lookup::cross::table_t& FormFactorManager::normalized_cross_table() noexcept {
    if (is_using_custom_form_factors()) {
        refresh_custom_state();
        return custom_tables->custom_normalized_cross_table;
    }
    return ::normalized_cross_table;
}

const lookup::atomic::table_t& FormFactorManager::raw_atomic_table() noexcept {
    if (is_using_custom_form_factors()) {
        refresh_custom_state();
        return custom_tables->custom_raw_atomic_table;
    }
    return ::raw_atomic_table;
}

const lookup::exv::table_t& FormFactorManager::raw_exv_table() noexcept {
    if (is_using_custom_form_factors()) {
        refresh_custom_state();
        return custom_tables->custom_raw_exv_table;
    }
    return ::raw_exv_table;
}

const lookup::cross::table_t& FormFactorManager::raw_cross_table() noexcept {
    if (is_using_custom_form_factors()) {
        refresh_custom_state();
        return custom_tables->custom_raw_cross_table;
    }
    return ::raw_cross_table;
}

void FormFactorManager::use_custom_form_factors(bool choice) {
    if (choice && !custom_tables) {
        std::vector<int> indices(settings::form_factor::max_ff_types);
        std::iota(indices.begin(), indices.end(), 0);
        set_custom_form_factors(indices);
    }
    _use_custom_form_factors = choice;
}

void FormFactorManager::set_custom_form_factors(std::vector<int> ff_indices) {
    assert(!ff_indices.empty() && "Custom form factors cannot be empty.");
    assert(ff_indices.size() <= form_factor::get_total_ff_count() && "Custom form factors cannot exceed the total number of available form factors.");

    _use_custom_form_factors = true;
    custom_tables = std::make_unique<_CustomTables>();
    custom_tables->active_count = ff_indices.size();
    std::copy(ff_indices.begin(), ff_indices.end(), custom_tables->ff_indices.begin());
    refresh();
}

void FormFactorManager::refresh() {
    _needs_refresh = true;
}

bool FormFactorManager::is_using_custom_form_factors() noexcept {
    return _use_custom_form_factors || !ExvTableManager::is_default();
}

void FormFactorManager::refresh_custom_state() {
    if (!_needs_refresh) {return;}
    assert(custom_tables && "Custom tables must be set before they can be used.");
    custom_tables->custom_raw_atomic_table = lookup::detail::generate_atomic_table<lookup::detail::RawFormFactorLookup>();
    custom_tables->custom_raw_cross_table = lookup::detail::generate_cross_table<lookup::detail::RawFormFactorLookup>();
    custom_tables->custom_raw_exv_table = lookup::detail::generate_exv_table();
    custom_tables->custom_normalized_atomic_table = lookup::detail::generate_atomic_table<lookup::detail::NormalizedFormFactorLookup>();
    custom_tables->custom_normalized_cross_table = lookup::detail::generate_cross_table<lookup::detail::NormalizedFormFactorLookup>();
    _needs_refresh = false;
}

void FormFactorManager::set_custom_form_factors(data::Molecule& molecule) {
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
    std::sort(ff_indices.begin(), ff_indices.end(), [&ff_counts](int a, int b) {
        return ff_counts[a] > ff_counts[b];
    });
    ff_indices.resize(settings::form_factor::max_ff_types);
    ff_indices.back() = static_cast<int>(form_factor::form_factor_t::OTHER); // OTHER will never be selected, so it is safe to assign it here

    {   // logging
        std::string log_msg = "Setting form factors based on detected molecular composition:";
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::EXCLUDED_VOLUME) + " (forced)";
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::WATER) + " (forced)";
        for (unsigned int i = 2; i < ff_indices.size()-1; ++i) {
            log_msg += "\n\t" + form_factor::to_string(static_cast<form_factor_t>(ff_indices[i])) + " with count " + std::to_string(ff_counts[ff_indices[i]]);
        }
        log_msg += "\n\t" + form_factor::to_string(form_factor::form_factor_t::OTHER) + " (forced)";
        logging::log(log_msg);
    }
    set_custom_form_factors(ff_indices);
}