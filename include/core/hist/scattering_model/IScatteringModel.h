// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/DistanceHistogram.h>
#include <math/Vector3.h>
#include <utility/Limit.h>

#include <optional>
#include <vector>

namespace ausaxs::hist {
    struct AmplitudeScatterers;

    /**
     * @brief Interface for a scattering model, which evaluates the distance histograms of a molecule together with
     *        its hydration model, and for IScatteringModelExv its excluded volume model, into the scattering intensity.
     *
     * A plain DistanceHistogram stores only the total pairwise distance distribution. A scattering
     * model instead keeps it split into three partial histograms, one per interaction type:
     *   - aa: atom-atom   — the solute with itself
     *   - aw: atom-water  — the solute with its hydration shell
     *   - ww: water-water — the hydration shell with itself
     *
     * The total histogram is simply their sum. Keeping the partials separate allows the
     * hydration-shell contribution to be rescaled independently of the solute (see
     * apply_water_scaling_factor()); this water-scaling factor is one of the parameters the fitter
     * optimizes. Each partial can also be transformed into its own scattering profile.
     */
    class IScatteringModel : public hist::DistanceHistogram {
        public:
            using hist::DistanceHistogram::DistanceHistogram;
            IScatteringModel() = default;
            IScatteringModel(const IScatteringModel&) = delete;
            IScatteringModel(IScatteringModel&&) noexcept = default;
            IScatteringModel& operator=(const IScatteringModel&) = delete;
            IScatteringModel& operator=(IScatteringModel&&) noexcept = default;
            ~IScatteringModel() override = default;

            /// @brief Get the partial distance histogram for atom-atom interactions.
            virtual const Distribution1D& get_aa_counts() const = 0;
            virtual Distribution1D& get_aa_counts() = 0; //< @copydoc get_aa_counts

            /// @brief Get the partial distance histogram for atom-water interactions.
            virtual const Distribution1D& get_aw_counts() const = 0;
            virtual Distribution1D& get_aw_counts() = 0; //< @copydoc get_aw_counts

            /// @brief Get the partial distance histogram for water-water interactions.
            virtual const Distribution1D& get_ww_counts() const = 0;
            virtual Distribution1D& get_ww_counts() = 0; //< @copydoc get_ww_counts

            /**
             * @brief Apply a scaling factor to the water partial distance histogram.
             */
            virtual void apply_water_scaling_factor(double k) = 0;

            /**
             * @brief Get the partial intensity profile for atom-atom interactions.
             */
            virtual ScatteringProfile get_profile_aa() const = 0;

            /**
             * @brief Get the partial intensity profile for atom-water interactions.
             */
            virtual ScatteringProfile get_profile_aw() const = 0;

            /**
             * @brief Get the partial intensity profile for water-water interactions.
             */
            virtual ScatteringProfile get_profile_ww() const = 0;

            /**
             * @brief Get the limits for the water scaling factor parameter. 
             *        This is intended to be used by the fitter to set correct limits. 
             */
            virtual Limit get_water_scaling_factor_limits() const {return {0, 10};}

            /**
             * @brief Determine if the structure represented by this histogram is highly ordered.
             */
            bool is_highly_ordered() const override;

            /**
             * @brief Evaluate the orientationally averaged scattering amplitude <A(q)> of @a scatterers at each of @a q, with the form factors and the 
             *        current parameters of this model. 
             *
             * This is the amplitude counterpart of the intensity I(q) = <|A(q)|^2>, and <A(q)>^2/I(q) is the beta(q) factor of the decoupling 
             * approximation, I(q) = I_0(q) [1 + beta(q) (S(q) - 1)]. Unlike the intensity, the amplitude depends on the point it is taken about, 
             * which must be the reference point of the particles in S(q).
             *
             * @param scatterers The scatterers of the molecule, from IHistogramManager::get_amplitude_scatterers() of the manager which built this model.
             * @param q The q values to evaluate at. The sums are over the binned distances, so any q value is supported.
             * @param centre The point the amplitude is taken about. Defaults to the centroid of the excess scattering length, see AmplitudeHistogram::scattering_centre.
             */
            std::vector<double> amplitude_transform(
                const AmplitudeScatterers& scatterers, const std::vector<double>& q, const std::optional<Vector3<double>>& centre = std::nullopt
            ) const;

        protected:
            /**
             * @brief The effective form factor g_c(q) of each class of the AmplitudeScatterers this model is built with: a scatterer of weight w in 
             *        class c adds w*g_c(q) to the amplitude. This is the amplitude each scatterer contributes to the intensity of this model.
             *        Throws by default, for the models without an amplitude.
             *
             * @param form_factors Whether the scatterers are resolved by form factor type, see AmplitudeScatterers.
             */
            virtual std::vector<double> amplitude_form_factors(bool form_factors, double q) const;
    };
}