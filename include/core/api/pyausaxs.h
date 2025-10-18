// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#ifdef WIN32
    #define API __declspec(dllexport)
#else
    #define API
#endif

extern "C" API void deallocate(int object_id, int* status);

extern "C" API int pdb_read(
    const char* filename,
    int* status
);

extern "C" API int pdb_get_data(
    int object_id,
    int** serial, const char*** name, const char*** altLoc, const char*** resName, const char** chainID, 
    int** resSeq, const char*** iCode, double** x, double** y, double** z, 
    double** occupancy, double** tempFactor, const char*** element, const char*** charge, 
    int* n_atoms, int* status
);

extern "C" API int data_read(
    const char* filename,
    int* status
);

extern "C" API int data_get_data(
    int object_id,
    double** q, double** I, double** Ierr, int* n_points,
    int* status
);

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
    double** ax, double** ay, double** az, double** aw, const char*** aform_factors, // atoms
    double** wx, double** wy, double** wz, double** ww,                              // waters
    int* an, int* wn, int* status
);