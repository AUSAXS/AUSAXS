// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/ScatteringAmplitudeNoExv.h>

namespace ausaxs::hist {
    /**
     * @brief The amplitude counterpart of ScatteringModelFraser. The form factors are those of ScatteringAmplitudeNoExv, less the excluded volume
     *        form factor of each atomic type scaled as in the intensity, crho cx(q) B_x(q) fx_t(q). The waters have no excluded volume.
     */
    class ScatteringAmplitudeFraser : public ScatteringAmplitudeNoExv {
        public:
            using ScatteringAmplitudeNoExv::ScatteringAmplitudeNoExv;

        protected:
            std::vector<double> form_factors(double q, const ScatteringParameters& params) const override;
    };
}
