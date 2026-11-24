// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/scattering_model/ScatteringParameters.h>
#include <math/Vector3.h>

#include <optional>
#include <vector>

namespace ausaxs::hist {
    /**
     * @brief The orientationally averaged scattering amplitude <A(q)> of a molecule and its hydration shell.
     *
     * The orientational average of the phase of a scatterer a distance r from the centre is sinc(qr), so
     *      <A(q)> = sum_k f_k(q) sum_{j in k} w_j sinc(q r_j)
     * where the scatterers are grouped by their form factor f_k, and w_j is the weight of scatterer j.
     */
    class ScatteringAmplitude {
        public:
            virtual ~ScatteringAmplitude() = default;

            /**
             * @brief Evaluate <A(q)> at each of @a q.
             *
             * @param q The q values to evaluate at.
             * @param centre The point the amplitude is taken about. Defaults to default_centre().
             * @param params The parameters of the scattering model.
             */
            std::vector<double> evaluate(
                const std::vector<double>& q, const std::optional<Vector3<double>>& centre = std::nullopt, const ScatteringParameters& params = {}
            ) const;

            /**
             * @brief The contrast-weighted centre of the molecule, i.e. the centroid of its scatterers weighted by their excess scattering length w_j f_k(0).
             *        The form factors are those of the scattering model, so the excluded volume is already subtracted from f_k(0) for the models that have one.
             *        This is the only centre where beta(q) has no q^2 term.
             */
            Vector3<double> default_centre(const ScatteringParameters& params = {}) const;

        protected:
            /**
             * @param groups The scatterers with their weights, grouped by form factor.
             */
            explicit ScatteringAmplitude(std::vector<detail::CompactCoordinates>&& groups);

            /**
             * @brief The excess form factor f_k(q) of each group k of scatterers with the model parameters @a params.
             */
            virtual std::vector<double> get_excess_form_factors(double q, const ScatteringParameters& params) const = 0;

        private:
            std::vector<detail::CompactCoordinates> groups;
    };
}
