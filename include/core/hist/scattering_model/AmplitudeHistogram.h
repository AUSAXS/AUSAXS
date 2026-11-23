// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <math/Vector3.h>

#include <functional>
#include <vector>

namespace ausaxs::hist {
    /**
     * @brief The scatterers of a molecule as a histogram manager weighs them, grouped into classes which share one effective form factor.
     *        A scatterer of weight w in class c adds w*g_c(q) to the amplitude, where the scattering model supplies g_c.
     *
     * The layout of the classes matches the scattering model the same manager builds:
     *  - without form factors: [atoms, waters], each scatterer with its own weight.
     *  - with form factors: one class per active form factor slot, with the waters in form_factor::water_bin, and every scatterer weighing 1.
     */
    struct AmplitudeScatterers {
        std::vector<detail::CompactCoordinates> classes;
        bool form_factors = false;

        static constexpr int atom_class = 0;  // the class of the atoms without form factors
        static constexpr int water_class = 1; // the class of the waters without form factors
    };

    /**
     * @brief The distances of the scatterers of each class from a centre, binned with the mean distance of each bin.
     *        This is the amplitude counterpart of the distance histogram: the orientational average of the phase of a scatterer a distance r from
     *        the centre is sinc(qr), so
     *              <A(q)> = sum_c g_c(q) sum_j w_j sinc(q r_j)
     *        Since each bin keeps the mean distance of its scatterers, the binning error is of second order in the bin width.
     */
    class AmplitudeHistogram {
        public:
            /**
             * @brief Bin the distances of @a scatterers from @a centre.
             */
            AmplitudeHistogram(const AmplitudeScatterers& scatterers, const Vector3<double>& centre);

            /**
             * @brief The centroid of the excess scattering length of @a scatterers, where @a g0 holds g_c(0) of each class.
             *        This is the default centre: it is the only one where beta(q) = <A(q)>^2/I(q) has no q^2 term.
             *        If the excess nearly cancels, as for a contrast-matched particle, it is ill-defined, and the absolute excess is used instead.
             */
            static Vector3<double> scattering_centre(const AmplitudeScatterers& scatterers, const std::vector<double>& g0);

            /**
             * @brief Evaluate <A(q)> at each of @a q, where @a form_factors(q) gives g_c(q) of each class.
             */
            std::vector<double> transform(const std::vector<double>& q, const std::function<std::vector<double>(double)>& form_factors) const;

        private:
            struct Class {std::vector<double> weights, distances;}; // only the occupied bins
            std::vector<Class> classes;
    };
}
