#pragma once

#include <hist/intensity_calculator/DistanceHistogram.h>
#include <utility/Limit.h>

namespace ausaxs::hist {
    class ICompositeDistanceHistogram : public hist::DistanceHistogram {
        public:
            using hist::DistanceHistogram::DistanceHistogram;
            ICompositeDistanceHistogram() = default;
            ICompositeDistanceHistogram(const ICompositeDistanceHistogram&) = default;
            ICompositeDistanceHistogram(ICompositeDistanceHistogram&&) noexcept = default;
            ICompositeDistanceHistogram& operator=(const ICompositeDistanceHistogram&) = default;
            ICompositeDistanceHistogram& operator=(ICompositeDistanceHistogram&&) noexcept = default;
            virtual ~ICompositeDistanceHistogram() = default;

            virtual const Distribution1D& get_aa_counts() const = 0;
            virtual Distribution1D& get_aa_counts() = 0;

            virtual const Distribution1D& get_aw_counts() const = 0;
            virtual Distribution1D& get_aw_counts() = 0;

            virtual const Distribution1D& get_ww_counts() const = 0;
            virtual Distribution1D& get_ww_counts() = 0;

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