// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <math/Vector3.h>

#include <vector>

namespace ausaxs::hist {
    /**
     * @brief Vector-Jacobian product of the raw Debye intensity with respect to the atomic coordinates.
     *
     * The raw intensity is the Debye sum of the atoms as weighted point scatterers,
     *      I(q) = \sum_{ij} w_i w_j sinc(q r_ij),
     * which the binned calculation approximates. Given the adjoint v(q) = dL/dI(q) at the values @a q, this returns
     *      dL/dr_i = \sum_q v(q) dI(q)/dr_i
     * for every atom, in the order of Molecule::iterate_atoms(). The derivative is that of the exact sum, not of its binned
     * approximation.
     *
     * The adjoint collapses the q axis before the pair loop: each pair only needs the scalar
     *      H(r) = 2 \sum_q v(q) (cos(qr) - sinc(qr))/r^2,
     * which is tabulated once per call, so the pair loop costs the same as a single distance histogram.
     * Hydrated molecules are not supported yet.
     */
    std::vector<Vector3<double>> debye_raw_vjp(const data::Molecule& molecule, const std::vector<double>& q, const std::vector<double>& v);
}
