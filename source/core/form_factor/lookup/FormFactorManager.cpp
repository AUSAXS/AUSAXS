#include <form_factor/lookup/FormFactorManager.h>
#include <container/Container2D.h>
#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/detail/FormFactorProductBase.h>
#include <form_factor/lookup/detail/LookupHelpers.h>

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
    CONST auto raw_exv_table           = lookup::detail::generate_exv_table();
    CONST auto raw_cross_table         = lookup::detail::generate_cross_table<lookup::detail::RawFormFactorLookup>();
    CONST auto normalized_atomic_table = lookup::detail::generate_atomic_table<lookup::detail::NormalizedFormFactorLookup>();
    CONST auto normalized_exv_table    = lookup::detail::generate_exv_table();
    CONST auto normalized_cross_table  = lookup::detail::generate_cross_table<lookup::detail::NormalizedFormFactorLookup>();
}

const lookup::atomic::table_t& FormFactorManager::normalized_atomic_table() noexcept {
    if (is_using_custom_form_factors()) {
        refresh_custom_state();
        return custom_tables->custom_normalized_atomic_table;
    }
    return ::normalized_atomic_table;
}

const lookup::exv::table_t& FormFactorManager::normalized_exv_table() noexcept {
    if (is_using_custom_form_factors()) {
        refresh_custom_state();
        return custom_tables->custom_normalized_exv_table;
    }
    return ::normalized_exv_table;
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
        std::vector<int> indices(form_factor::get_count_without_excluded_volume());
        std::iota(indices.begin(), indices.end(), 0);
        set_custom_form_factors(indices);
    }
    _use_custom_form_factors = choice;
}

void FormFactorManager::set_custom_form_factors(std::vector<int> ff_indices) {
    assert(!ff_indices.empty() && "Custom form factors cannot be empty.");
    assert(ff_indices.size() <= form_factor::get_count_without_excluded_volume() && "Custom form factors cannot exceed the total number of available form factors.");

    _use_custom_form_factors = true;
    custom_tables = std::make_unique<_CustomTables>();
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
    custom_tables->custom_normalized_exv_table = lookup::detail::generate_exv_table();
    _needs_refresh = false;
}