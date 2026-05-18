// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/DistanceHistogram.h>
#include <utility/Limit.h>

namespace ausaxs::hist {
    /**
     * @brief Interface for a composite distance histogram.
     *
     * A plain DistanceHistogram stores only the total pairwise distance distribution. A *composite*
     * histogram instead keeps it split into three partial histograms, one per interaction type:
     *   - aa: atom-atom   — the solute with itself
     *   - aw: atom-water  — the solute with its hydration shell
     *   - ww: water-water — the hydration shell with itself
     *
     * The total histogram is simply their sum. Keeping the partials separate allows the
     * hydration-shell contribution to be rescaled independently of the solute (see
     * apply_water_scaling_factor()); this water-scaling factor is one of the parameters the fitter
     * optimizes. Each partial can also be transformed into its own scattering profile.
     */
    class ICompositeDistanceHistogram : public hist::DistanceHistogram {
        public:
            using hist::DistanceHistogram::DistanceHistogram;
            ICompositeDistanceHistogram() = default;
            ICompositeDistanceHistogram(const ICompositeDistanceHistogram&) = default;
            ICompositeDistanceHistogram(ICompositeDistanceHistogram&&) noexcept = default;
            ICompositeDistanceHistogram& operator=(const ICompositeDistanceHistogram&) = default;
            ICompositeDistanceHistogram& operator=(ICompositeDistanceHistogram&&) noexcept = default;
            virtual ~ICompositeDistanceHistogram() = default;

            /// @brief Get the partial distance histogram for atom-atom interactions.
            virtual const Distribution1D& get_aa_counts() const = 0;
            virtual Distribution1D& get_aa_counts() = 0; //< @copydoc get_aa_counts

            /// @brief Get the partial distance histogram for atom-water interactions.
            virtual const Distribution1D& get_aw_counts() const = 0;
            virtual Distribution1D& get_aw_counts() = 0; //< @copydoc get_aw_counts

            /// @brief Get the partial distance histogram for water-water interactions.
            virtual const Distribution1D& get_ww_counts() const = 0;
            virtual Distribution1D& get_ww_counts() = 0; //< @copydoc get_ww_counts

            /**
             * @brief Apply a scaling factor to the water partial distance histogram.
             */
            virtual void apply_water_scaling_factor(double k) = 0;

            /**
             * @brief Get the partial intensity profile for atom-atom interactions.
             */
            virtual ScatteringProfile get_profile_aa() const = 0;

            /**
             * @brief Get the partial intensity profile for atom-water interactions.
             */
            virtual ScatteringProfile get_profile_aw() const = 0;

            /**
             * @brief Get the partial intensity profile for water-water interactions.
             */
            virtual ScatteringProfile get_profile_ww() const = 0;

            /**
             * @brief Get the limits for the water scaling factor parameter. 
             *        This is intended to be used by the fitter to set correct limits. 
             */
            virtual Limit get_water_scaling_factor_limits() const {return {0, 10};}
    };
}