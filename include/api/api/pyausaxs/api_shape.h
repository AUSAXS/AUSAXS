// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

/**
 * @brief The scattering intensity of a homogeneous body made of occupied cubic lattice cells.
 *        This uses a lattice transform to compute the Debye sum in O(N log N) time, where N is the number of cells.
 *
 * @param x, y, z The cell centres. They must all lie on one lattice of the given spacing, and no two may coincide.
 * @param n_cells The number of cells.
 * @param spacing The lattice spacing in Å.
 * @param q The q values to evaluate, in Å^-1.
 * @param I Output: the intensity at each q.
 * @param n_q The number of q values.
 */
extern "C" API void shape_debye_userq(
    const double* x, const double* y, const double* z, int n_cells, double spacing,
    const double* q, double* I, int n_q,
    int* status
);