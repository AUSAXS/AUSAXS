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
 * The three Debye families below differ in what they include, and share their q axis:
 *
 * - molecule_debye: the molecule's configured model (histogram manager and excluded volume), including its hydration
 *   shell if it has one. With the default Simple model every atom carries the one Gaussian form factor of the Debye
 *   transform, so I(q) keeps an overall exp(-q^2); the other exv models apply their per-species form factors instead.
 * - molecule_debye_raw: the plain binned Debye sum of the atoms and any waters as point scatterers with their weights:
 *   no excluded volume, and the exp(-q^2) divided out.
 * - molecule_debye_exact: as _raw, but summed over all atom pairs without distance binning, and over the atoms only:
 *   a hydration shell is ignored. Slow; meant as a reference for _raw on a molecule without waters.
 *
 * The variants without a q argument evaluate on the default q axis from settings::axes::qmin to settings::axes::qmax and
 * return those q values; the _userq variants evaluate at the given q values. The returned id owns the arrays and is
 * released with deallocate.
 */
extern "C" API int molecule_debye(
    int molecule_id,
    double** q, double** I, int* n_points, 
    int* status
);

extern "C" API void molecule_debye_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
);

extern "C" API int molecule_debye_raw(
    int molecule_id,
    double** q, double** I, int* n_points,
    int* status
);

extern "C" API void molecule_debye_raw_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
);

extern "C" API int molecule_debye_exact(
    int molecule_id,
    double** q, double** I, int* n_points,
    int* status
);

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