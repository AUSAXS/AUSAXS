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
 * @brief Evaluate the I(q) for a given set of q-values and a set of coordinates.
 *        The coordinates are assumed to be in the form of a list of x, y, z, and w values.
 *        The return status will be non-zero if an error occurred.
 *        The resulting I(q) values will be stored in the given array.
 */
extern "C" API void evaluate_sans_debye(double* _q, double* _x, double* _y, double* _z, double* _w, int _nq, int _nc, int* _return_status, double* _return_Iq);