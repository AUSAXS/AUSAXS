// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorManager.h>
#include <hist/scattering_model/ScatteringModelFraserBase.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of ScatteringModelSimple, this class also uses individual form factors for each atomic type.
     *        Unique exluded volume form factors are used for each atomic type.
     *        For more information, see ScatteringModelSimple.
     */
    class ScatteringModelFraser : public ScatteringModelFraserBase<form_factor::lookup::table_t, form_factor::lookup::table_t, form_factor::lookup::table_t> {
        public: 
            using ScatteringModelFraserBase::ScatteringModelFraserBase;

            const form_factor::lookup::table_t& get_ff_table() const override {
                return form_factor::manager::get_active_product_tables()->raw_atomic_table;
            }

            const form_factor::lookup::table_t& get_ffax_table() const override {
                return form_factor::manager::get_active_product_tables()->raw_cross_table;
            }

            const form_factor::lookup::table_t& get_ffxx_table() const override {
                return form_factor::manager::get_active_product_tables()->raw_exv_table;
            }

        protected:
            /**
             * @brief The atomic form factors less the excluded volume form factor of each type, scaled as in the intensity:
             *        g_t(q) = B_a(q)*f_t(q) - rho*cx(q)*B_x(q)*fx_t(q), and c_w*f_w(q) for the waters, which carry no excluded volume.
             */
            std::vector<double> amplitude_form_factors(bool form_factors, double q) const override;
    };
    static_assert(supports_nothrow_move_v<ScatteringModelFraser>, "ScatteringModelFraser should be nothrow move constructible");
}