// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/ExvTable.h>
#include <form_factor/ExvFormFactor.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::form_factor {
    class ExvTableManager {
        public:
            static observer_ptr<const constants::exv::detail::ExvSet> get_current_exv_table();
            static constexpr detail::ExvFormFactorSet get_current_exv_form_factor_set();
            static void set_custom_exv_table(const constants::exv::detail::ExvSet& set);
            static bool is_default();

        private:
            static inline bool _use_custom_exv_table = false;
            static inline std::unique_ptr<constants::exv::detail::ExvSet> custom_exv_tables;
            static const detail::ExvFormFactorSet& _nonconstexpr_get_current_exv_form_factor_set();
    };
}

constexpr ausaxs::form_factor::detail::ExvFormFactorSet ausaxs::form_factor::ExvTableManager::get_current_exv_form_factor_set() {
    if (std::is_constant_evaluated()) {
        return detail::ExvFormFactorSet(constants::exv::standard);
    }
    return _nonconstexpr_get_current_exv_form_factor_set();
}