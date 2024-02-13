#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <container/Container1D.h>
#include <container/Container2D.h>
#include <container/Container3D.h>

namespace hist {
    class ICompositeDistanceHistogramExv : public ICompositeDistanceHistogram {
        public:
            ICompositeDistanceHistogramExv() = default;
            using ICompositeDistanceHistogram::ICompositeDistanceHistogram;
            virtual ~ICompositeDistanceHistogramExv() = default;

            /**
             * @brief Apply a scaling factor to the excluded volume partial distance histogram.
             */
            virtual void apply_excluded_volume_scaling_factor(double k) = 0;

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual ScatteringProfile get_profile_ax() const = 0;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual ScatteringProfile get_profile_xx() const = 0;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual ScatteringProfile get_profile_wx() const = 0;
    };
}