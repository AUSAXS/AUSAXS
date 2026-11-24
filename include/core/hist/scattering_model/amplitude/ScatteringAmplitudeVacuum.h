// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/amplitude/ScatteringAmplitude.h>

namespace ausaxs::hist {
    /**
     * @brief The amplitude in vacuum, i.e. without an excluded volume; the counterpart of ScatteringModelNoExv. The atoms of each form factor type t have the form factor B_a(q) f_t(q), with
     *        the atomic Debye-Waller factor B_a, and the waters have cw f_w(q).
     */
    class ScatteringAmplitudeVacuum : public ScatteringAmplitude {
        public:
            /**
             * @param groups The scatterers of each active form factor type with unit weights, with the waters in form_factor::water_bin.
             */
            explicit ScatteringAmplitudeVacuum(std::vector<detail::CompactCoordinates>&& groups);

        protected:
            container::Container2D<double> get_excess_form_factors(const std::vector<double>& q, const ScatteringParameters& params) const override;
    };
}
