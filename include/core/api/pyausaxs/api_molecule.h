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
    double** ax, double** ay, double** az, double** aw, const char*** aff,
    double** wx, double** wy, double** wz, double** ww,
    int* an, int* wn, int* status
);

extern "C" API void molecule_hydrate(
    int molecule_id,
    int* status
);

extern "C" API int molecule_distance_histogram(
    int molecule_id,
    double** aa, double** aw, double** ww,
    int* n_bins, double* delta_r, int* status
);

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

extern "C" int molecule_debye_raw(
    int molecule_id,
    double** q, double** I, int* n_points,
    int* status
);

extern "C" void molecule_debye_raw_userq(
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