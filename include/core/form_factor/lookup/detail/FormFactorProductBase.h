// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <settings/FormFactorSettings.h>
#include <constants/Constants.h>

#include <array>

namespace ausaxs::form_factor::lookup::detail {
    /**
     * @brief Generate an atomic form factor product table.
     * @tparam FormFactorLookup A type providing a static `get(form_factor_t)` method.
     */
    template<typename FormFactorLookup>
    const form_factor::lookup::atomic::table_t generate_atomic_table(const std::array<int, settings::form_factor::max_ff_types>& ff_indices) {
        form_factor::lookup::atomic::table_t table;
        for (unsigned int i = 0; i < ff_indices.size(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                table.index(i, j) = FormFactorProduct(
                    FormFactorLookup::get(static_cast<form_factor_t>(ff_indices[i])), 
                    FormFactorLookup::get(static_cast<form_factor_t>(ff_indices[j]))
                );
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(
                FormFactorLookup::get(static_cast<form_factor_t>(ff_indices[i])), 
                FormFactorLookup::get(static_cast<form_factor_t>(ff_indices[i]))
            );
        }

        return table;
    }

    /**
     * @brief Generate an excluded volume form factor product table (exv-exv).
     *        This is a symmetric table.
     * @param use_default_table If true, always use the default EXV set.
     */
    const inline form_factor::lookup::exv::table_t generate_exv_table(const std::array<int, settings::form_factor::max_ff_types>& ff_indices, bool use_default_table = false) {
        auto exv_set = use_default_table ? ExvTableManager::get_default_exv_form_factor_set() : ExvTableManager::get_current_exv_form_factor_set();

        form_factor::lookup::exv::table_t table;
        for (unsigned int i = start_index_for_explicit_exv(); i < ff_indices.size(); ++i) {
            for (unsigned int j = start_index_for_explicit_exv(); j < i; ++j) {
                table.index(i, j) = FormFactorProduct(
                    exv_set.get(static_cast<form_factor_t>(ff_indices[i])), 
                    exv_set.get(static_cast<form_factor_t>(ff_indices[j]))
                );
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(
                exv_set.get(static_cast<form_factor_t>(ff_indices[i])), 
                exv_set.get(static_cast<form_factor_t>(ff_indices[i]))
            );
        }
        return table;
    }

    /**
     * @brief Generate a cross form factor product table (atomic-exv).
     * @tparam AtomicFormFactorLookup A type providing a static `get(form_factor_t)` method for atomic form factors.
     * @param use_default_table If true, always use the default EXV set.
     */
    template<typename AtomicFormFactorLookup>
    const form_factor::lookup::cross::table_t generate_cross_table(const std::array<int, settings::form_factor::max_ff_types>& ff_indices, bool use_default_table = false) {
        auto exv_set = use_default_table ? ExvTableManager::get_default_exv_form_factor_set() : ExvTableManager::get_current_exv_form_factor_set();

        form_factor::lookup::cross::table_t table;
        for (unsigned int i = 0; i < ff_indices.size(); ++i) {
            for (unsigned int j = start_index_for_explicit_exv(); j < ff_indices.size(); ++j) {
                table.index(i, j) = FormFactorProduct(
                    AtomicFormFactorLookup::get(static_cast<form_factor_t>(ff_indices[i])), 
                    exv_set.get(static_cast<form_factor_t>(ff_indices[j]))
                );
            }
        }
        return table;
    }
}