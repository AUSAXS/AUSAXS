// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/DistributionFwd.h>

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        Unique exluded volume form factors are used for each atomic type.
     *        For more information, see CompositeDistanceHistogram.
     */
    template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
    class CompositeDistanceHistogramFFExplicitBase : public CompositeDistanceHistogramFFAvgBase<AAFormFactorTableType> {
        public: 
            CompositeDistanceHistogramFFExplicitBase();
            CompositeDistanceHistogramFFExplicitBase(const CompositeDistanceHistogramFFExplicitBase&);
            CompositeDistanceHistogramFFExplicitBase(CompositeDistanceHistogramFFExplicitBase&&) noexcept;
            CompositeDistanceHistogramFFExplicitBase& operator=(CompositeDistanceHistogramFFExplicitBase&&) noexcept;
            CompositeDistanceHistogramFFExplicitBase& operator=(const CompositeDistanceHistogramFFExplicitBase&);
            virtual ~CompositeDistanceHistogramFFExplicitBase() override;

            /**
             * @brief Create a new unweighted composite distance histogram with form factors.
             *        The same distance histogram is used for aa, ax, and xx interactions (with different form factor tables).
             *        Similarly, the same histogram is used for aw and wx interactions.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions (also used for ax and xx).
             * @param p_aw The partial distance histogram for atom-water interactions (also used for wx).
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             */
            CompositeDistanceHistogramFFExplicitBase(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot
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
             */
            CompositeDistanceHistogramFFExplicitBase(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot
            );

            const AAFormFactorTableType get_ffaa_table() const;
            virtual const AXFormFactorTableType& get_ffax_table() const = 0;
            virtual const XXFormFactorTableType& get_ffxx_table() const = 0;

            /**
             * @brief Get the excluded volume scaling factor.
             *
             * @param cx The scaling factor for the excluded volume.
             * @param q The scattering vector.
             */
            static double exv_factor(double q, double cx);

        protected:
            // @copydoc CompositeDistanceHistogramFFAvgBase::exv_factor(double) const
            double exv_factor(double q) const override;

            //#################################//
            //###           CACHE           ###//
            //#################################//

            mutable struct {
                // cached sinqd vals for each form factor combination
                // indexing as [ff1][ff2] for 3D, [ff1] for 2D
                // Note: aa, ax, xx share the same underlying histogram, so we only cache sinqd_aa
                // Similarly, aw and wx share the same histogram
                mutable struct {
                    container::Container3D<double> aa;
                    container::Container2D<double> aw;
                    container::Container1D<double> ww;
                    bool valid = false;
                } sinqd;
            } exv_cache;

        private:
            void cache_refresh_intensity_profiles(bool sinqd_changed, bool cw_changed, bool cx_changed) const override;
            void cache_refresh_sinqd() const override;
    };
}