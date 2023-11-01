#pragma once

#include <hist/intensity_calculator/DistanceHistogram.h>

namespace hist {
    struct ICompositeDistanceHistogram : public hist::DistanceHistogram {
        using hist::DistanceHistogram::DistanceHistogram;
        virtual ~ICompositeDistanceHistogram() = default;

        virtual const std::vector<double>& get_aa_counts() const = 0;
        virtual std::vector<double>& get_aa_counts() = 0;

        virtual const std::vector<double>& get_aw_counts() const = 0;
        virtual std::vector<double>& get_aw_counts() = 0;

        virtual const std::vector<double>& get_ww_counts() const = 0;
        virtual std::vector<double>& get_ww_counts() = 0;

        virtual void apply_water_scaling_factor(double k) = 0;

        virtual const ScatteringProfile get_profile_aa() const = 0;

        virtual const ScatteringProfile get_profile_aw() const = 0;

        virtual const ScatteringProfile get_profile_ww() const = 0;
    };
}