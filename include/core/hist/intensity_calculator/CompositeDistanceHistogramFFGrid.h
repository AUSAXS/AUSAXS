// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridBase.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        The excluded volume is approximated using a space-filling grid of spheres, filling the volume of the molecule. 
     *        This approach adds a substantial overhead to the calculations, but should give a more accurate representation of the excluded volume.
     */
    class CompositeDistanceHistogramFFGrid : public CompositeDistanceHistogramFFGridBase {
        // needs to move the cached state out of a separately-evaluated FFGrid instance
        friend class CompositeDistanceHistogramFFGridScalableExv;
        public:
            CompositeDistanceHistogramFFGrid(CompositeDistanceHistogramFFGrid&&) noexcept;
            CompositeDistanceHistogramFFGrid& operator=(CompositeDistanceHistogramFFGrid&&) noexcept;
            ~CompositeDistanceHistogramFFGrid() override;

            /**
             * @brief Create a weighted grid-based composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot_aa The total distance histogram for everything except the grid. This is only used to extract the bin centers.
             * @param p_tot_ax The total distance histogram for the cross terms. This is only used to extract the bin centers. 
             * @param p_tot_xx The total distance histogram for the grid only. Calculations involving this grid must use unique bin centers due to the highly ordered grid structure. 
             */
            CompositeDistanceHistogramFFGrid(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot_aa,
                hist::WeightedDistribution1D&& p_tot_ax,
                hist::WeightedDistribution1D&& p_tot_xx
            );

        protected:
            double exv_factor(double q) const override;

            /**
             * @brief The grid excluded volume is independent data, so it is read straight out of the distance
             *        histograms - but using the grid's own sinc(qd) tables rather than the atomic one.
             */
            void cache_refresh_sinqd_exv() const override;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFGrid>, "CompositeDistanceHistogramFFGrid should support nothrow move semantics.");
}
