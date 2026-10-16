// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API int molecule_from_file(
    const char* filename,
    int* status
);

extern "C" API int molecule_from_pdb_id(
    int pdb_id,
    int* status
);

extern "C" API int molecule_from_arrays(
    double* x, double* y, double* z, double* w, int n_atoms,
    int* status
);

extern "C" API int molecule_get_data(
    int molecule_id,
    double** ax_out, double** ay_out, double** az_out, double** aw_out, const char*** aform_factors_out,
    double** wx_out, double** wy_out, double** wz_out, double** ww_out,
    int* na, int* nw, int* status
);

extern "C" API void molecule_hydrate(
    int molecule_id,
    int* status
);

extern "C" API int molecule_distance_histogram(
    int molecule_id,
    double** aa, double** aw, double** ww, double** axis, int* n_bins,
    int* status
);

/**
 * The total distance histogram transformed by molecule_debye_raw: the atoms and any waters as point scatterers with their weights, without any excluded 
 * volume. The counts sum the weight products over all ordered pairs, so the self-terms are included in the first bin.
 */
extern "C" API int molecule_distance_histogram_raw(
    int molecule_id,
    double** counts, double** axis, int* n_bins,
    int* status
);

/**
 * The Debye intensity of the molecule's configured model, including its hydration shell if it has one.
 */
extern "C" API int molecule_debye(
    int molecule_id,
    double** q, double** I, int* n_points, 
    int* status
);

/**
 * As molecule_debye, but evaluated at the given q values.
 */
extern "C" API void molecule_debye_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
);

/**
 * The plain binned Debye sum of the atoms and any waters as point scatterers with their weights.
 */
extern "C" API int molecule_debye_raw(
    int molecule_id,
    double** q, double** I, int* n_points,
    int* status
);

/**
 * As molecule_debye_raw, but evaluated at the given q values.
 */
extern "C" API void molecule_debye_raw_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
);

/**
 * As molecule_debye_raw, but summed over all pairs without distance binning. Exists only as a reference for molecule_debye_raw; far too slow for 
 * practical use.
 */
extern "C" API int molecule_debye_exact(
    int molecule_id,
    double** q, double** I, int* n_points,
    int* status
);

/**
 * As molecule_debye_exact, but evaluated at the given q values.
 */
extern "C" API void molecule_debye_exact_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
);

extern "C" API int molecule_debye_fit(
    int molecule_id, int data_id,
    int* status
);

extern "C" API void molecule_clear_hydration(
    int molecule_id,
    int* status
);

extern "C" API void molecule_Rg(
    int molecule_id,
    double* Rg,
    int* status
);