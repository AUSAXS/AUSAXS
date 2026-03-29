// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/observer_ptr.h>
#include <container/ArrayContainer2D.h>
#include <form_factor/lookup/FormFactorLookupFwd.h>

#include <memory>

namespace ausaxs::form_factor {
    class FormFactorManager {
        struct _CustomTables{
            std::array<int, form_factor::get_count_without_excluded_volume()> ff_indices;
            lookup::exv::table_t    custom_raw_exv_table;
            lookup::cross::table_t  custom_raw_cross_table;
            lookup::atomic::table_t custom_raw_atomic_table;
            lookup::exv::table_t    custom_normalized_exv_table;
            lookup::cross::table_t  custom_normalized_cross_table;
            lookup::atomic::table_t custom_normalized_atomic_table;
        };

        public:
            static observer_ptr<const _CustomTables> get_custom_tables() noexcept;
            static constexpr std::array<int, form_factor::get_count_without_excluded_volume()> get_ff_indices() noexcept;
            static const lookup::atomic::table_t& normalized_atomic_table() noexcept;
            static const lookup::exv::table_t& normalized_exv_table() noexcept;
            static const lookup::cross::table_t& normalized_cross_table() noexcept;
            static const lookup::atomic::table_t& raw_atomic_table() noexcept;
            static const lookup::exv::table_t& raw_exv_table() noexcept;
            static const lookup::cross::table_t& raw_cross_table() noexcept;
            static void use_custom_form_factors(bool choice);
            static void set_custom_form_factors(std::vector<int> ff_indices);
            static void refresh();

        private:
            static inline bool _needs_refresh = false;
            static inline bool _use_custom_form_factors = false;
            static inline std::unique_ptr<_CustomTables> custom_tables;
            static bool is_using_custom_form_factors() noexcept;

            static void refresh_custom_state();
    };
}

constexpr std::array<int, ausaxs::form_factor::get_count_without_excluded_volume()> ausaxs::form_factor::FormFactorManager::get_ff_indices() noexcept {
    auto generator = []() {
        std::array<int, form_factor::get_count_without_excluded_volume()> indices{};
        for (unsigned int i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }
        return indices;
    };

    if (std::is_constant_evaluated()) {return generator();}
    if (_use_custom_form_factors) {return custom_tables->ff_indices;}
    return generator();
}