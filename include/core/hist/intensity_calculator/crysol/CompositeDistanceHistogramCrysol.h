// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <utility/TypeTraits.h>
#include <utility/observer_ptr.h>

namespace ausaxs::hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the CRYSOL excluded volume fitting.
     *        The form factor tables are shared with CompositeDistanceHistogramFFExplicit, and only the excluded volume scaling G(q)
     *        and its fitting limits differ. Like CRYSOL, this always uses the Traube volumes: constructing it sets settings::exv::exv_set accordingly.
     */
    class CompositeDistanceHistogramCrysol : public CompositeDistanceHistogramFFExplicit {
        public:
            CompositeDistanceHistogramCrysol() = default;

            /**
             * @brief Create a new unweighted composite distance histogram with form factors.
             *        The same distance histogram is used for aa, ax, and xx interactions (with different form factor tables).
             *        Similarly, the same histogram is used for aw and wx interactions.
             *
             * @param p_aa The partial distance histogram for atom-atom interactions (also used for ax and xx).
             * @param p_aw The partial distance histogram for atom-water interactions (also used for wx).
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             * @param molecule The molecule the histograms were calculated from. Its average displaced volume per atom determines G(q).
             */
            CompositeDistanceHistogramCrysol(
                hist::Distribution3D&& p_aa,
                hist::Distribution2D&& p_aw,
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot,
                observer_ptr<const data::Molecule> molecule
            );

            /**
             * @brief Create a new weighted composite distance histogram with form factors.
             *        The same distance histogram is used for aa, ax, and xx interactions (with different form factor tables).
             *        Similarly, the same histogram is used for aw and wx interactions.
             *
             * @param p_aa The partial distance histogram for atom-atom interactions (also used for ax and xx).
             * @param p_aw The partial distance histogram for atom-water interactions (also used for wx).
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used to extract the bin centers.
             * @param molecule The molecule the histograms were calculated from. Its average displaced volume per atom determines G(q).
             */
            CompositeDistanceHistogramCrysol(
                hist::Distribution3D&& p_aa,
                hist::Distribution2D&& p_aw,
                hist::Distribution1D&& p_ww,
                hist::WeightedDistribution1D&& p_tot,
                observer_ptr<const data::Molecule> molecule
            );

            Limit get_excluded_volume_scaling_factor_limits() const override;

            /**
             * @brief Get the excluded volume scaling factor.
             *
             * @param cx The scaling factor for the excluded volume.
             * @param q The scattering vector.
             * @param avg_displaced_V The average displaced volume per atom.
             */
            static double exv_factor(double q, double cx, double avg_displaced_V);

            double average_displaced_V = 0;

        protected:
            double exv_factor(double q) const override;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramCrysol>, "CompositeDistanceHistogramCrysol should be nothrow move constructible");
}
