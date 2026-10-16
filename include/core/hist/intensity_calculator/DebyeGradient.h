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
     * which the binned calculation approximates. Given the loss gradient v(q) = dL/dI(q) at the values @a q, this returns
     *      dL/dr_i = \sum_q v(q) dI(q)/dr_i
     * for every atom. The derivative is that of the exact sum, not of its binned approximation.
     */
    std::vector<Vector3<double>> debye_raw_vjp(const data::Molecule& molecule, const std::vector<double>& q, const std::vector<double>& v);
}
