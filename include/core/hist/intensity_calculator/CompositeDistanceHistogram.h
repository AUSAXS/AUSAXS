#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/distribution/DistributionFwd.h>
#include <hist/distribution/Distribution1D.h>
#include <container/Container1D.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions.
     *        This allows scaling the water contribution to the total intensity profile for use in e.g. fitting. 
     */
    class CompositeDistanceHistogram : public ICompositeDistanceHistogram {
        public: 
            CompositeDistanceHistogram(const CompositeDistanceHistogram&);
            CompositeDistanceHistogram(CompositeDistanceHistogram&&) noexcept;
            CompositeDistanceHistogram& operator=(CompositeDistanceHistogram&&) noexcept;
            virtual ~CompositeDistanceHistogram() override;

            /**
             * @brief Create an unweighted composite distance histogram.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram.
             */
            CompositeDistanceHistogram(
                hist::Distribution1D&& p_aa, 
                hist::Distribution1D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::Distribution1D&& p_tot
            );

            CompositeDistanceHistogram(
                hist::Distribution1D&& p_aa, 
                hist::Distribution1D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot
            );

            virtual const hist::Distribution1D& get_aa_counts() const override;
            virtual hist::Distribution1D& get_aa_counts() override;
            virtual const Distribution1D& get_aw_counts() const override;
            virtual Distribution1D& get_aw_counts() override;
            virtual const Distribution1D& get_ww_counts() const override;
            virtual Distribution1D& get_ww_counts() override;

            virtual void apply_water_scaling_factor(double k) override;
            virtual ScatteringProfile get_profile_aa() const override;
            virtual ScatteringProfile get_profile_aw() const override;
            virtual ScatteringProfile get_profile_ww() const override;

        private:
            mutable struct {Distribution1D aa, aw, ww;} distance_profiles;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogram>, "CompositeDistanceHistogram should support nothrow move semantics.");
}