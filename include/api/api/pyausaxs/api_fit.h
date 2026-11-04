// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API int fit_get_fit_info(
    int fit_id,
    const char*** pars, double** pvals, double** perr_min, double** perr_max, int* n_pars,
    double* chi_squared, int* dof,
    int* status
);

extern "C" API int fit_get_fit_curves(
    int fit_id,
    double** q, double** I_data, double** I_err, double** I_model, int* n_points,
    int* status
);
/**
 * @brief Get the radius of gyration of the fitted model: the q -> 0 limit of its scattering curve, with the excluded
 *        volume subtracted and the hydration shell at its fitted contrast. NaN if the fit did not record one.
 */
extern "C" API void fit_get_Rg(
    int fit_id,
    double* Rg,
    int* status
);
