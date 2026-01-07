// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

#if CONSTEXPR_TABLES
    #define FF_CONST constexpr
#else
    #define FF_CONST const
#endif

namespace ausaxs::form_factor::lookup::detail {
    /**
     * @brief Generate an atomic form factor product table.
     * 
     * @tparam FormFactorLookup A type providing a static `get(form_factor_t)` method.
     */
    template<typename FormFactorLookup>
    FF_CONST form_factor::lookup::atomic::table_t generate_atomic_table() {
        form_factor::lookup::atomic::table_t table;
        for (unsigned int i = 0; i < form_factor::get_count(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                table.index(i, j) = FormFactorProduct(
                    FormFactorLookup::get(static_cast<form_factor_t>(i)), 
                    FormFactorLookup::get(static_cast<form_factor_t>(j))
                );
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(
                FormFactorLookup::get(static_cast<form_factor_t>(i)), 
                FormFactorLookup::get(static_cast<form_factor_t>(i))
            );
        }
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
    FF_CONST TableType generate_exv_table(const form_factor::detail::ExvFormFactorSet& set) {
        TableType table;
        for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                table.index(i, j) = FormFactorProduct(
                    set.get(static_cast<form_factor_t>(i)), 
                    set.get(static_cast<form_factor_t>(j))
                );
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(
                set.get(static_cast<form_factor_t>(i)), 
                set.get(static_cast<form_factor_t>(i))
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
    FF_CONST TableType generate_cross_table(const form_factor::detail::ExvFormFactorSet& exv_set) {
        TableType table;
        for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < form_factor::get_count_without_excluded_volume(); ++j) {
                table.index(i, j) = FormFactorProduct(
                    AtomicFormFactorLookup::get(static_cast<form_factor_t>(i)), 
                    exv_set.get(static_cast<form_factor_t>(j))
                );
            }
        }
        return table;
    }

    /**
     * @brief Get the appropriate exv table based on current settings.
     * 
     * @tparam TableType The table type.
     * @tparam TableGenerator A callable that takes an ExvFormFactorSet and returns the table.
     * @param default_table The default table to return when using default settings.
     * @param custom_table The custom table to return when using custom settings.
     * @param generator The table generator function.
     */
    template<typename TableType, typename TableGenerator>
    const TableType& get_exv_table_for_settings(
        const TableType& default_table,
        const TableType& custom_table,
        TableGenerator&& generator
    ) {
        if (settings::molecule::exv_set == settings::molecule::ExvSet::Default) {
            return default_table;
        }
        switch (settings::molecule::exv_set) {
            case settings::molecule::ExvSet::Traube: {
                static auto table = generator(form_factor::detail::ExvFormFactorSet(constants::exv::Traube));
                return table;
            }
            case settings::molecule::ExvSet::Voronoi_explicit_H: {
                static auto table = generator(form_factor::detail::ExvFormFactorSet(constants::exv::Voronoi_explicit_H));
                return table;
            }
            case settings::molecule::ExvSet::Voronoi_implicit_H: {
                static auto table = generator(form_factor::detail::ExvFormFactorSet(constants::exv::Voronoi_implicit_H));
                return table;
            }
            case settings::molecule::ExvSet::MinimumFluctutation_explicit_H: {
                static auto table = generator(form_factor::detail::ExvFormFactorSet(constants::exv::MinimumFluctuation_explicit_H));
                return table;
            }
            case settings::molecule::ExvSet::MinimumFluctutation_implicit_H: {
                static auto table = generator(form_factor::detail::ExvFormFactorSet(constants::exv::MinimumFluctuation_implicit_H));
                return table;
            }
            case settings::molecule::ExvSet::vdw: {
                static auto table = generator(form_factor::detail::ExvFormFactorSet(constants::exv::vdw));
                return table;
            }
            case settings::molecule::ExvSet::Custom: {
                return custom_table;
            }
            default: 
                return default_table;
        }
    }
}