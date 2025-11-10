// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicitBase.h>
#include <hist/intensity_calculator/foxs/FormFactorFoXS.h>
#include <form_factor/xray/PrecalculatedExvFormFactorProduct.h>

namespace ausaxs::hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that uses the same form factors as FoXS. 
     */
    class CompositeDistanceHistogramFoXS : public CompositeDistanceHistogramFFExplicitBase<form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t> {
        public: 
            using CompositeDistanceHistogramFFExplicitBase::CompositeDistanceHistogramFFExplicitBase;
            ~CompositeDistanceHistogramFoXS() override = default;

            CompositeDistanceHistogramFoXS(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot
            );

            CompositeDistanceHistogramFoXS(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot
            );

            Limit get_excluded_volume_scaling_factor_limits() const override;

            const form_factor::storage::atomic::table_t& get_ff_table() const override {
                return ff_aa_table;
            }

            const form_factor::storage::cross::table_t& get_ffax_table() const override {
                return ff_ax_table;
            }

            const form_factor::storage::exv::table_t& get_ffxx_table() const override {
                return ff_xx_table;
            }

            /**
             * @brief Get the excluded volume scaling factor.
             *
             * @param cx The scaling factor for the excluded volume.
             * @param q The scattering vector.
             */
            static double exv_factor(double q, double cx);

        protected:
            double exv_factor(double q) const override;
            form_factor::storage::atomic::table_t ff_aa_table = form_factor::foxs::storage::atomic::generate_table();
            form_factor::storage::cross::table_t ff_ax_table  = form_factor::foxs::storage::cross::generate_table();
            form_factor::storage::exv::table_t ff_xx_table    = form_factor::foxs::storage::exv::generate_table();
    };
}