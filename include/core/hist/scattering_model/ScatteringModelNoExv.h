// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorLookupFwd.h>
#include <hist/scattering_model/ScatteringModelExvBase.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A scattering model where all scatterers are taken to be in vacuum, i.e. there is no excluded volume contribution. 
     */
    class ScatteringModelNoExv : public ScatteringModelExvBase<form_factor::lookup::table_t> {
        public:
            using ScatteringModelExvBase::ScatteringModelExvBase;

        protected:
            const form_factor::lookup::table_t& get_ff_table() const override;

            void cache_refresh_sinqd_exv() const override {}
            void cache_refresh_intensity_exv(const std::vector<double>& /* void */, bool /* void */, bool /* void */) const override {}

            /**
             * @brief Evaluate the intensity directly at arbitrary q values, without using precomputed tables.
             */
            std::vector<double> debye_transform_direct(const std::vector<double>& q) const override;
    };
    static_assert(supports_nothrow_move_v<ScatteringModelNoExv>, "ScatteringModelNoExv should support nothrow move semantics.");
}
