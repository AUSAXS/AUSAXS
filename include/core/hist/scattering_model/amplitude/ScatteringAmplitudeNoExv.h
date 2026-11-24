// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/amplitude/ScatteringAmplitude.h>

namespace ausaxs::hist {
    /**
     * @brief The amplitude counterpart of ScatteringModelNoExv. The atoms of each form factor type t have the form factor B_a(q) f_t(q), with
     *        the atomic Debye-Waller factor B_a, and the waters have cw f_w(q).
     */
    class ScatteringAmplitudeNoExv : public ScatteringAmplitude {
        public:
            /**
             * @param groups The scatterers of each active form factor type with unit weights, with the waters in form_factor::water_bin.
             */
            explicit ScatteringAmplitudeNoExv(std::vector<detail::CompactCoordinates>&& groups);

        protected:
            std::vector<double> form_factors(double q, const ScatteringParameters& params) const override;
    };
}
