#pragma once

#include <hist/intensity_calculator/DistanceHistogram.h>

namespace hist {
    class ICompositeDistanceHistogram : public hist::DistanceHistogram {
        public:
            using hist::DistanceHistogram::DistanceHistogram;
            virtual ~ICompositeDistanceHistogram() = default;

            virtual const Distribution1D& get_aa_counts() const = 0;
            virtual Distribution1D& get_aa_counts() = 0;

            virtual const Distribution1D& get_aw_counts() const = 0;
            virtual Distribution1D& get_aw_counts() = 0;

            virtual const Distribution1D& get_ww_counts() const = 0;
            virtual Distribution1D& get_ww_counts() = 0;

            virtual void apply_water_scaling_factor(double k) = 0;
            void reset_water_scaling_factor() {apply_water_scaling_factor(1);}

            virtual const ScatteringProfile get_profile_aa() const = 0;

            virtual const ScatteringProfile get_profile_aw() const = 0;

            virtual const ScatteringProfile get_profile_ww() const = 0;
    };
}