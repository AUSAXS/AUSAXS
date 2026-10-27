// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

/**
 * @brief Test the integration of the C++ code with the Python code.
 *        This will simply increment the value of the given pointer by 1.
 */
extern "C" API void test_integration(int* test_value);

/**
 * @brief Evaluate the I(q) for a given set of q-values and a set of coordinates.
 *        The coordinates are assumed to be in the form of a list of x, y, z, and w values.
 *        The return status will be non-zero if an error occurred.
 *        The resulting I(q) values will be stored in the given array.
 *        No form factors, excluded volume or implicit hydrogens are used; the global settings are left as they were.
 */
extern "C" API void debye_no_ff(
    double* _q, double* _x, double* _y, double* _z, double* _w, 
    int _nq, int _nc, double* _return_Iq, int* _return_status
);
