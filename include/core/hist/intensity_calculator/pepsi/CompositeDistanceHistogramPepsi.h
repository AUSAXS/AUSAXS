// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace ausaxs::hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the Pepsi-SAXS excluded volume fitting. 
     */
    class CompositeDistanceHistogramPepsi : public CompositeDistanceHistogramCrysol {
        public:
            using CompositeDistanceHistogramCrysol::CompositeDistanceHistogramCrysol;
            ~CompositeDistanceHistogramPepsi() override = default;

            Limit get_excluded_volume_scaling_factor_limits() const override;

            /**
             * @brief Get the excluded volume scaling factor.
             *
             * @param cx The scaling factor for the excluded volume.
             * @param q The scattering vector.
             */
            static double exv_factor(double q, double cx);

        protected:
            double exv_factor(double q) const override;
    };
}