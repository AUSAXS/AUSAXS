// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API void ff_get_five_gaussian_coefficients(
    const char* element, 
    double* a, double* b, double* c,
    int* status
);

extern "C" API void ff_get_current_exv_volume(
    const char* element,
    double* volume,
    int* status
);

extern "C" API int ff_valid_form_factor_types(
    const char*** types,
    int* n_types,
    int* status
);
