// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

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

extern "C" API int pdb_debye_fit(
    int pdb_id, int data_id,
    int* status
);

/**
 * @brief Fit the named symmetry to the chains of a PDB and return the symmetry-decomposed structure: the reference chain plus the fitted 
 *        symmetry copies, each atom tagged with its copy index (0 = reference chain, 1..N = generated copies). `rmsd` receives the residual 
 *        RMSD of the fit (Å), a direct measure of how well the assembly obeys the claimed symmetry.
 */
extern "C" API int pdb_decompose_symmetry(
    int pdb_id, const char* symmetry,
    double** x, double** y, double** z, int** copy_index, int* n_atoms,
    double* rmsd, int* status
);