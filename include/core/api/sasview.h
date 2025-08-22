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
 * @brief Start an iterative fit.
 *
 * This will initialize the fit process, converting the given data and PDB coordinates into a
 * dataset and a molecule, respectively. The fit will be performed in steps, allowing SasView to control the fitting process.
*/
void iterative_fit_start(
    double* data_q, double* data_I, double* data_Ierr, int n_data,
    double* pdb_x,  double* pdb_y,  double* pdb_z, 
    const char** atom_names, const char** residue_names, const char** elements, 
    int n_pdb, int* return_status
);

/**
 * @brief Calculate the scattering intensity for the current parameters.
 */
void iterative_fit_step(double* pars, double* return_I, int* return_status);

/**
 * @brief Finish the iterative fit process.
 *        This will return the final intensity, and write out the model to disk. 
 */
void iterative_fit_finish(double* pars, double* return_I, int* return_status);

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
