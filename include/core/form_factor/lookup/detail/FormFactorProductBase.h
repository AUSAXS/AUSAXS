// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

#if CONSTEXPR_TABLES
    #define CONST constexpr
#else
    #define CONST const
#endif
namespace ausaxs::form_factor::lookup::detail {
    /**
     * @brief Generate an atomic form factor product table.
     * 
     * @tparam FormFactorLookup A type providing a static `get(form_factor_t)` method.
     */
    template<typename FormFactorLookup>
    CONST form_factor::lookup::atomic::table_t generate_atomic_table() {
        auto ff_indices = FormFactorManager::get_ff_indices();

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

        // add excluded volume ff for the exv models using the same form factor for all
        int i = form_factor::get_count()-1;
        for (unsigned int j = 0; j < ff_indices.size(); ++j) {
            table.index(i, j) = FormFactorProduct(
                FormFactorLookup::get(static_cast<form_factor_t>(ff_indices[j])), 
                FormFactorLookup::get(form_factor_t::EXCLUDED_VOLUME)
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = FormFactorProduct(
            FormFactorLookup::get(form_factor_t::EXCLUDED_VOLUME), 
            FormFactorLookup::get(form_factor_t::EXCLUDED_VOLUME)
        );

        return table;
    }

    /**
     * @brief Generate an excluded volume form factor product table (exv-exv).
     *        This is a symmetric table.
     * 
     * @tparam TableType The table type to generate.
     * @param set The excluded volume form factor set to use.
     */
    template<typename TableType>
    CONST TableType generate_exv_table() {
        auto ff_indices = FormFactorManager::get_ff_indices();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        TableType table;
        for (unsigned int i = 0; i < ff_indices.size(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
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
     * 
     * @tparam TableType The table type to generate.
     * @tparam AtomicFormFactorLookup A type providing a static `get(form_factor_t)` method for atomic form factors.
     * @param exv_set The excluded volume form factor set to use.
     */
    template<typename TableType, typename AtomicFormFactorLookup>
    CONST TableType generate_cross_table() {
        auto ff_indices = FormFactorManager::get_ff_indices();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        TableType table;
        for (unsigned int i = 0; i < ff_indices.size(); ++i) {
            for (unsigned int j = 0; j < ff_indices.size(); ++j) {
                table.index(i, j) = FormFactorProduct(
                    AtomicFormFactorLookup::get(static_cast<form_factor_t>(ff_indices[i])), 
                    exv_set.get(static_cast<form_factor_t>(ff_indices[j]))
                );
            }
        }
        return table;
    }
}