// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief An instantiation of the CompositeDistanceHistogramFFAvgBase class that uses the default precalculated form factor table.
     *        A single averaged excluded volume form factor is used for all atoms.
     *        For more information, see CompositeDistanceHistogramFFAvgBase.
     *
     * Every excluded volume "atom" of this model sits on top of a real atom and carries the same averaged
     * form factor, so the excluded volume distance data is not independent: it is a fixed linear combination
     * of the atomic distance data, scaled by the average displaced charge. It is therefore never stored;
     * see cache_refresh_sinqd_exv.
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
             * @brief Derive the excluded volume sinqd values from the atomic ones.
             *
             * With one averaged excluded volume "atom" per real atom, and Z the average displaced charge:
             *     ax(a) = Z * sum_b [aa(a,b) + aa(b,a)],   xx = Z^2 * sum_ab aa(a,b),   wx = 2Z * sum_a aw(a)
             * so no excluded volume distance histogram is needed - the sums above are exact, not approximate.
             * ax and wx already count both directions, per the convention of the sinqd cache.
             */
            void cache_refresh_sinqd_exv() const override;

            double Z_exv_avg = 0; // the average excluded volume charge displaced by a single atom
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFAvg>, "CompositeDistanceHistogramAvg should support nothrow move semantics.");
}
