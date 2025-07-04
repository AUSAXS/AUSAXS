// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#ifdef WIN32
    #define API __declspec(dllexport)
#else
    #define API
#endif

/**
 * @brief Test the integration of the C++ code with the Python code.
 *        This will simply increment the value of the given pointer by 1.
 */
extern "C" API void test_integration(int* test_value);

/**
 * @brief Fit the scattering intensity for a structure to the given data.
 *        The _pdb_type is the form factor type of the atoms in the structure; see FormFactorType.h for possible values.
 *        The return status will be non-zero if an error occurred.
 *        The resulting I(q) values will be stored in the given array.
 */
extern "C" API void fit_saxs_future(
    double* _data_q, double* _data_I, double* _data_Ierr, int _n_data,
    double* _pdb_x,  double* _pdb_y,  double* _pdb_z,     int _pdb_type, int _n_pdb,
    double* _return_I, int* _return_status
);

/**
 * @brief Fit the scattering intensity for a structure to the given data.
 *        The return status will be non-zero if an error occurred.
 *        The resulting I(q) values will be stored in the given array.
 */
extern "C" API void fit_saxs(
    double* _data_q, double* _data_I, double* _data_Ierr, int _n_data,
    double* _pdb_x,  double* _pdb_y,  double* _pdb_z, const char** atom_names, 
    const char** residue_names, const char** elements, int _n_pdb,
    double* _return_I, int* _return_status
);

/**
 * @brief Evaluate the I(q) for a given set of q-values and a set of coordinates.
 *        The coordinates are assumed to be in the form of a list of x, y, z, and w values.
 *        The return status will be non-zero if an error occurred.
 *        The resulting I(q) values will be stored in the given array.
 */
extern "C" API void evaluate_sans_debye(
    double* _q, double* _x, double* _y, double* _z, double* _w, 
    int _nq, int _nc, double* _return_Iq, int* _return_status
);
