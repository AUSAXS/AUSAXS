// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

/**
 * @brief The scattering intensity of a homogeneous body made of occupied cubic lattice cells.
 *        Uses the lattice transform and a binned Debye sum.
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

/**
 * @brief The scattering intensity of a homogeneous body on an orthogonal rectangular lattice.
 *
 * Each cell has volume spacing_x*spacing_y*spacing_z and unit contrast. Like shape_debye_userq,
 * the cell form factor is a volume-equivalent isotropic Gaussian, not a rectangular-box form factor.
 * Pair distances are binned and the Debye result may be interpolated at the requested q values.
 *
 * @param x, y, z The unique cell centres, on a lattice with the given axis spacings.
 * @param n_cells The number of cells.
 * @param spacing_x, spacing_y, spacing_z Positive lattice spacings in Å.
 * @param q The q values to evaluate, in Å^-1.
 * @param I Output: the intensity at each q.
 * @param n_q The number of q values.
 */
extern "C" API void shape_debye_rectangular_userq(
    const double* x, const double* y, const double* z, int n_cells,
    double spacing_x, double spacing_y, double spacing_z,
    const double* q, double* I, int n_q,
    int* status
);