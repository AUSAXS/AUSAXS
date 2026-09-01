// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <utility/TypeTraits.h>

#include <memory>

namespace ausaxs::hist {
    /**
     * @brief An instantiation of the CompositeDistanceHistogramFFAvgBase class that uses the default precalculated form factor table.
     */
    class CompositeDistanceHistogramFFAvg : public CompositeDistanceHistogramFFAvgBase<form_factor::lookup::table_t> {
        using CompositeDistanceHistogramFFAvgBase::CompositeDistanceHistogramFFAvgBase;

        public:
            /**
             * @brief Synthesize the averaged excluded volume distributions into @a p_aa and @a p_aw, and construct a histogram from the result.
             *
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             * @param Z_exv_avg The average excluded volume charge displaced by a single atom.
             */
            [[nodiscard]] static std::unique_ptr<CompositeDistanceHistogramFFAvg> with_averaged_exv(
                hist::Distribution3D&& p_aa,
                hist::Distribution2D&& p_aw,
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot,
                double Z_exv_avg
            );

            /**
             * @brief Weighted overload of with_averaged_exv. 
             *        @a p_tot is only used to extract the bin centers.
             */
            [[nodiscard]] static std::unique_ptr<CompositeDistanceHistogramFFAvg> with_averaged_exv(
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
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFAvg>, "CompositeDistanceHistogramAvg should support nothrow move semantics.");
}
