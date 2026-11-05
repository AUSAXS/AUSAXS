// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/ScatteringModelExvBase.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief The atoms and waters with their atomic form factors, without any excluded volume. This is ExvMethod::None for the form factor-resolved
     *        managers, see settings::internal_state::prefer_explicit_form_factors. The excluded volume profiles ax, xx and wx are always zero.
     */
    class ScatteringModelNoExv : public ScatteringModelExvBase<form_factor::lookup::table_t> {
        public:
            using ScatteringModelExvBase::ScatteringModelExvBase;

        protected:
            const form_factor::lookup::table_t& get_ff_table() const override {
                return form_factor::manager::get_active_product_tables()->raw_atomic_table;
            }

            void cache_refresh_sinqd_exv() const override {}
            void cache_refresh_intensity_exv(const std::vector<double>& /* void */, bool /* void */, bool /* void */) const override {}

            // without an excluded volume, the model holds at any q
            std::vector<double> debye_transform_beyond_axis(const std::vector<double>& q) const override {return atomic_debye_transform(q);}
    };
    static_assert(supports_nothrow_move_v<ScatteringModelNoExv>, "ScatteringModelNoExv should support nothrow move semantics.");
}
