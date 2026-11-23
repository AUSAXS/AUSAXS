// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>

#include <vector>

namespace ausaxs::hist {
    /**
     * @brief Evaluate the orientationally averaged scattering amplitude <A(q)> of the molecule, including its hydration shell.
     *
     * Each scatterer j has an effective form factor g_j(q), the amplitude it adds to the excess scattering. Averaged over all orientations,
     * the phase of a scatterer at distance r_j from the centre reduces to sinc(q r_j), so
     *      <A(q)> = sum_j g_j(q) sinc(q r_j)
     * This is the amplitude counterpart of the Debye intensity I(q) = <|A(q)|^2>, and <A(q)>^2 / I(q) is the beta(q) factor of the
     * Kotlarchyk-Chen decoupling approximation, I(q) = I_0(q) [1 + beta(q) (S(q) - 1)].
     *
     * Unlike the intensity, the amplitude depends on the choice of centre beyond the Guinier region. The centre used here is the centroid of the
     * excess scattering length, weighting each scatterer by g_j(0), which is the natural centre of the scattering particle.
     *
     * The effective form factors follow settings::exv::exv_method, with the conventions of the intensity of the plain histogram managers at their
     * default, unfitted parameters. The sums are exact, so any q value is supported.
     * Only the Simple, None, Vacuum and Fraser excluded volume models are supported, and molecules with symmetries are not; both throw.
     */
    std::vector<double> average_amplitude(const data::Molecule& molecule, const std::vector<double>& q_vals);
}
