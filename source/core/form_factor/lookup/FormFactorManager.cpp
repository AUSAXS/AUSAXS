#pragma once

#include <form_factor/lookup/FormFactorManager.h>
#include <container/Container2D.h>
#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/detail/FormFactorProductBase.h>
#include <form_factor/lookup/detail/LookupHelpers.h>

#include <numeric>

using namespace ausaxs;
using namespace ausaxs::form_factor; 

observer_ptr<const FormFactorManager::_CustomTables> FormFactorManager::get_all_tables() noexcept {
    return custom_tables.get();
}

namespace {
    // Default tables. These _must_ be initialized statically to ensure default settings. 
    constexpr auto raw_atomic_table        = lookup::detail::generate_atomic_table<lookup::detail::RawFormFactorLookup>();
    constexpr auto raw_exv_table           = lookup::detail::generate_exv_table<lookup::exv::table_t>();
    constexpr auto raw_cross_table         = lookup::detail::generate_cross_table<lookup::exv::table_t, lookup::detail::RawFormFactorLookup>();
    constexpr auto normalized_atomic_table = lookup::detail::generate_atomic_table<lookup::detail::NormalizedFormFactorLookup>();
    constexpr auto normalized_exv_table    = lookup::detail::generate_exv_table<lookup::exv::table_t>();
    constexpr auto normalized_cross_table  = lookup::detail::generate_cross_table<lookup::exv::table_t, lookup::detail::NormalizedFormFactorLookup>();
}
const lookup::atomic::table_t& FormFactorManager::normalized_atomic_table() const noexcept {
    if (_use_custom_form_factors) {
        refresh_custom_state();
        return custom_tables->custom_normalized_atomic_table;
    }
    return ::normalized_atomic_table;
}

const lookup::exv::table_t& FormFactorManager::normalized_exv_table() const noexcept {
    if (_use_custom_form_factors) {
        refresh_custom_state();
        return custom_tables->custom_normalized_exv_table;
    }
    return ::normalized_exv_table;
}

const lookup::cross::table_t& FormFactorManager::normalized_cross_table() const noexcept {
    if (_use_custom_form_factors) {
        refresh_custom_state();
        return custom_tables->custom_normalized_cross_table;
    }
    return ::normalized_cross_table;
}

const lookup::atomic::table_t& FormFactorManager::raw_atomic_table() const noexcept {
    if (_use_custom_form_factors) {
        refresh_custom_state();
        return custom_tables->custom_raw_atomic_table;
    }
    return ::raw_atomic_table;
}

const lookup::exv::table_t& FormFactorManager::raw_exv_table() const noexcept {
    if (_use_custom_form_factors) {
        refresh_custom_state();
        return custom_tables->custom_raw_exv_table;
    }
    return ::raw_exv_table;
}

const lookup::cross::table_t& FormFactorManager::raw_cross_table() const noexcept {
    if (_use_custom_form_factors) {
        refresh_custom_state();
        return custom_tables->custom_raw_cross_table;
    }
    return ::raw_cross_table;
}

void FormFactorManager::use_custom_form_factors(bool choice) {
    if (!custom_tables) {throw std::logic_error("Custom form factors must be set before they can be used.");}
    _use_custom_form_factors = choice;
}

void FormFactorManager::set_custom_form_factors(std::vector<int> ff_indices) {
    custom_tables = std::make_unique<_CustomTables>();
    _needs_refresh = true;
    assert(!custom_form_factors.empty() && "Custom form factors cannot be empty.");
    assert(custom_form_factors.size() <= form_factor::get_count() && "Custom form factors cannot exceed the total number of available form factors.");
}

const std::array<int, form_factor::get_count_without_excluded_volume()>& FormFactorManager::get_ff_indices() noexcept {
    if (_use_custom_form_factors) {
        assert(custom_tables && "Custom tables must be set before they can be used.");
        return custom_tables->ff_indices;
    }

    static std::array<int, form_factor::get_count_without_excluded_volume()> indices;
    std::iota(indices.begin(), indices.end(), 0);
    return indices;
}

void FormFactorManager::refresh_custom_state() {
    if (!_needs_refresh) {return;}
    assert(custom_tables && "Custom tables must be set before they can be used.");
    custom_tables->custom_raw_atomic_table = lookup::detail::generate_atomic_table<lookup::detail::RawFormFactorLookup>();
    custom_tables->custom_raw_cross_table = lookup::detail::generate_cross_table<lookup::cross::table_t, lookup::detail::RawFormFactorLookup>();
    custom_tables->custom_raw_exv_table = lookup::detail::generate_exv_table<lookup::exv::table_t>();
    custom_tables->custom_normalized_atomic_table = lookup::detail::generate_atomic_table<lookup::detail::NormalizedFormFactorLookup>();
    custom_tables->custom_normalized_cross_table = lookup::detail::generate_cross_table<lookup::cross::table_t, lookup::detail::NormalizedFormFactorLookup>();
    custom_tables->custom_normalized_exv_table = lookup::detail::generate_exv_table<lookup::exv::table_t>();
    _needs_refresh = false;
}