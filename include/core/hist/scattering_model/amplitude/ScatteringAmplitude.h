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
     *        This is the amplitude counterpart of the intensity I(q) = <|A(q)|^2> of a scattering model, and is created by the same histogram
     *        manager, see IHistogramManager::get_scattering_amplitude(). Each subclass is the counterpart of one scattering model.
     *
     * The orientational average of the phase of a scatterer a distance r from the centre is sinc(qr), so
     *      <A(q)> = sum_k f_k(q) sum_{j in k} w_j sinc(q r_j)
     * where the scatterers are grouped by their form factor f_k, and w_j is the weight of scatterer j.
     * <A(q)>^2/I(q) is the beta(q) factor of the decoupling approximation, I(q) = I_0(q) [1 + beta(q) (S(q) - 1)].
     * Unlike the intensity, the amplitude depends on the point it is taken about, which must be the reference point of the particles in S(q).
     *
     * The distances from the centre are binned with the mean distance of each bin, so the binning error is of second order in the bin width,
     * and any q value is supported. Evaluating the amplitude costs O(N), and never requires the distance histograms.
     */
    class ScatteringAmplitude {
        public:
            virtual ~ScatteringAmplitude() = default;

            /**
             * @brief Evaluate <A(q)> at each of @a q.
             *
             * @param q The q values to evaluate at.
             * @param centre The point the amplitude is taken about. Defaults to default_centre().
             * @param params The parameters of the scattering model, for example IScatteringModel::get_parameters() of a fitted model.
             */
            std::vector<double> evaluate(
                const std::vector<double>& q, const std::optional<Vector3<double>>& centre = std::nullopt, const ScatteringParameters& params = {}
            ) const;

            /**
             * @brief The centroid of the scattering length of the molecule, i.e. of its scatterers weighted by w_j f_k(0).
             *        This is the only centre where beta(q) has no q^2 term.
             *        If the scattering length nearly cancels, as for a contrast-matched particle, the centroid of its absolute value is used instead.
             */
            Vector3<double> default_centre(const ScatteringParameters& params = {}) const;

        protected:
            /**
             * @param groups The scatterers with their weights, grouped by form factor.
             */
            explicit ScatteringAmplitude(std::vector<detail::CompactCoordinates>&& groups);

            /**
             * @brief The form factor f_k(q) of each group k of scatterers with the model parameters @a params.
             */
            virtual std::vector<double> form_factors(double q, const ScatteringParameters& params) const = 0;

        private:
            std::vector<detail::CompactCoordinates> groups;
    };
}
