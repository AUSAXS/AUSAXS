#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the Pepsi-SAXS excluded volume fitting. 
     */
    class CompositeDistanceHistogramPepsi : public CompositeDistanceHistogramCrysol {
        public:
            using CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol;
            ~CompositeDistanceHistogramPepsi() override = default;

            Limit get_excluded_volume_scaling_factor_limits() const override;

        protected:
            double exv_factor(double q) const override;
    };
}