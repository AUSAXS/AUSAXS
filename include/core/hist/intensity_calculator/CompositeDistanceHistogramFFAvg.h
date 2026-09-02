// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief The average excluded volume model, where all dummy exv atoms are placed on top of the real atoms and share a single average form factor.
     */
    class CompositeDistanceHistogramFFAvg : public CompositeDistanceHistogramFFAvgBase<form_factor::lookup::table_t> {
        public:
            CompositeDistanceHistogramFFAvg() = default;

            /**
             * @brief Create an unweighted composite distance histogram with an averaged excluded volume.
             *
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             * @param Z_exv_avg The average excluded volume charge displaced by a single atom.
             */
            CompositeDistanceHistogramFFAvg(
                hist::Distribution3D&& p_aa,
                hist::Distribution2D&& p_aw,
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot,
                double Z_exv_avg
            );

            /**
             * @brief Create a weighted composite distance histogram with an averaged excluded volume.
             *        @a p_tot is only used to extract the bin centers.
             */
            CompositeDistanceHistogramFFAvg(
                hist::Distribution3D&& p_aa,
                hist::Distribution2D&& p_aw,
                hist::Distribution1D&& p_ww,
                hist::WeightedDistribution1D&& p_tot,
                double Z_exv_avg
            );

        protected:
            const form_factor::lookup::table_t& get_ff_table() const override {
                return form_factor::manager::get_active_product_tables()->raw_atomic_table;
            }

            /**
             * @brief No-op: the dummy atoms sit on top of the real atoms, so the excluded volume terms share the
             *        atomic sinqd cache. cache_refresh_intensity_exv contracts it with the exv form factors directly.
             */
            void cache_refresh_sinqd_exv() const override {}

            void cache_refresh_intensity_exv(const std::vector<double>& cx, bool cw_changed, bool cx_changed) const override;

            double Z_exv_avg = 0; // the average excluded volume charge displaced by a single atom
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFAvg>, "CompositeDistanceHistogramAvg should support nothrow move semantics.");
}
