// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/amplitude/ScatteringAmplitudeVacuum.h>

namespace ausaxs::hist {
    /**
     * @brief The amplitude counterpart of ScatteringModelFraser. The form factors are those of ScatteringAmplitudeVacuum, less the excluded volume
     *        form factor of each atomic type scaled as in the intensity, crho cx(q) B_x(q) fx_t(q). The waters have no excluded volume.
     */
    class ScatteringAmplitudeFraser : public ScatteringAmplitudeVacuum {
        public:
            using ScatteringAmplitudeVacuum::ScatteringAmplitudeVacuum;

        protected:
            std::vector<double> get_excess_form_factors(double q, const ScatteringParameters& params) const override;
    };
}
