// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/amplitude/ScatteringAmplitude.h>

namespace ausaxs::hist {
    /**
     * @brief The amplitude counterpart of ScatteringModelSimple. Every scatterer is a point of its weight with the form factor exp(-q^2/2),
     *        since the intensity is exp(-q^2) sum_ij w_i w_j sinc(q r_ij). The waters are scaled by cw.
     */
    class ScatteringAmplitudeSimple : public ScatteringAmplitude {
        public:
            /**
             * @param atoms The atoms with their weights, including any excluded volume of the Simple model.
             * @param waters The waters with their weights.
             */
            ScatteringAmplitudeSimple(detail::CompactCoordinates&& atoms, detail::CompactCoordinates&& waters);

        protected:
            container::Container2D<double> get_excess_form_factors(const std::vector<double>& q, const ScatteringParameters& params) const override;
    };
}
