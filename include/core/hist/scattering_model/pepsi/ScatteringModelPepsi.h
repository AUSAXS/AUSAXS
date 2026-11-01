// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/scattering_model/ScatteringModelFraser.h>
#include <hist/scattering_model/crysol/ScatteringModelCrysol.h>

namespace ausaxs::hist {
    /**
     * @brief An alternative to ScatteringModelFraser that mimics the Pepsi-SAXS excluded volume fitting. 
     */
    class ScatteringModelPepsi : public ScatteringModelCrysol {
        public:
            using ScatteringModelCrysol::ScatteringModelCrysol;
            ~ScatteringModelPepsi() override = default;

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