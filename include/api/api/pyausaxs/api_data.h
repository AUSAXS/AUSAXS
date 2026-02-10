// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API int data_read(
    const char* filename,
    int* status
);

extern "C" API int data_get_data(
    int object_id,
    double** q, double** I, double** Ierr, int* n_points,
    int* status
);