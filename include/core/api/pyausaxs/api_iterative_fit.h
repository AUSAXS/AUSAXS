// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API int iterative_fit_init(
    int molecule_id, 
    int* status
);

extern "C" API int iterative_fit_init_userq(
    int molecule_id, 
    double* q, int n_points,
    int* status
);

extern "C" API void iterative_fit_evaluate(
    int iterative_fit_id, 
    double* pars, int n_pars, 
    double** return_I, int* n_points,
    int* status
);

extern "C" API void iterative_fit_evaluate_userq(
    int iterative_fit_id, 
    double* pars, int n_pars, 
    double* q, double* I, int n_points,
    int* status
);