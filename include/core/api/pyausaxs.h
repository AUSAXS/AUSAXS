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