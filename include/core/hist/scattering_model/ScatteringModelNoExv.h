// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorLookupFwd.h>
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
            const form_factor::lookup::table_t& get_ff_table() const override;

            void cache_refresh_sinqd_exv() const override {}
            void cache_refresh_intensity_exv(const std::vector<double>& /* void */, bool /* void */, bool /* void */) const override {}

            /**
             * @brief Without an excluded volume, the model holds at any q: the atomic and hydration terms aa + aw + ww, with the form factors
             *        evaluated at each q value rather than taken from the tables over the default q axis.
             *        The water scaling factor and the atomic Debye-Waller factor are applied.
             */
            std::vector<double> debye_transform_beyond_axis(const std::vector<double>& q) const override;
    };
    static_assert(supports_nothrow_move_v<ScatteringModelNoExv>, "ScatteringModelNoExv should support nothrow move semantics.");
}
