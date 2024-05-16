#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the CRYSOL excluded volume fitting. 
     */
    class CompositeDistanceHistogramCrysol : public CompositeDistanceHistogramFFExplicit {
        public:
            using CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit;
            ~CompositeDistanceHistogramCrysol() override = default;

            Limit get_excluded_volume_scaling_factor_limits() const override;

        protected:
            double exv_factor(double q) const override;
    };
}