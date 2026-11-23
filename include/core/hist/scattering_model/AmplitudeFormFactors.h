// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/ScatteringParameters.h>
#include <math/Vector3.h>

#include <memory>
#include <optional>
#include <vector>

namespace ausaxs::hist {
    struct AmplitudeScatterers;

    /**
     * @brief The form factors a scattering model gives its scatterers, separate from its distance histograms.
     *        This is all the orientationally averaged amplitude needs from the model, so it can be evaluated without calculating the histograms.
     *        It holds no parameters of its own: they are passed in, either the current ones of a model or the defaults.
     */
    class AmplitudeFormFactors {
        public:
            virtual ~AmplitudeFormFactors() = default;

            /**
             * @brief The effective form factor g_c(q) of each class of AmplitudeScatterers with the parameters @a params: a scatterer of weight w in
             *        class c adds w*g_c(q) to the amplitude. This is the amplitude each scatterer contributes to the intensity of the model.
             */
            virtual std::vector<double> evaluate(const ScatteringParameters& params, double q) const = 0;

            /**
             * @brief Whether the classes are resolved by form factor type, see AmplitudeScatterers.
             */
            virtual bool form_factors() const = 0;
    };

    /**
     * @brief The form factors of the scattering model of settings::exv::exv_method, or nullptr if it has no amplitude.
     *        This is the only mapping from the excluded volume models to their form factors; the models built by the histogram managers carry
     *        the same ones (see IScatteringModel::amplitude_transform).
     */
    std::shared_ptr<const AmplitudeFormFactors> make_amplitude_form_factors();

    /**
     * @brief Evaluate the orientationally averaged scattering amplitude <A(q)> of @a scatterers at each of @a q.
     *
     * This is the amplitude counterpart of the intensity I(q) = <|A(q)|^2>, and <A(q)>^2/I(q) is the beta(q) factor of the decoupling
     * approximation, I(q) = I_0(q) [1 + beta(q) (S(q) - 1)]. Unlike the intensity, the amplitude depends on the point it is taken about,
     * which must be the reference point of the particles in S(q).
     *
     * @param form_factors The form factors of the model. Throws if this is nullptr, i.e. the model has no amplitude.
     * @param scatterers The scatterers of the molecule, from IHistogramManager::get_amplitude_scatterers().
     * @param q The q values to evaluate at. The sums are over the binned distances, so any q value is supported.
     * @param params The parameters of the model.
     * @param centre The point the amplitude is taken about. Defaults to the centroid of the excess scattering length, see AmplitudeHistogram::scattering_centre.
     */
    std::vector<double> amplitude_transform(
        const AmplitudeFormFactors* form_factors, const AmplitudeScatterers& scatterers, const std::vector<double>& q,
        const ScatteringParameters& params, const std::optional<Vector3<double>>& centre
    );
}
