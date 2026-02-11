// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/pyausaxs/api_data.h>
#include <api/pyausaxs/api_fit.h>
#include <api/pyausaxs/api_iterative_fit.h>
#include <api/pyausaxs/api_molecule.h>
#include <api/pyausaxs/api_pdb.h>
#include <api/pyausaxs/api_rigidbody.h>
#include <api/pyausaxs/api_settings.h>

// extern "C" API int map_read(
//     const char* filename,
//     int* status
// );

// extern "C" API void map_get_slice(
//     int map_id,
//     double z_position,
//     double** slice_data,
//     int* width, int* height,
//     int* status
// );

// extern "C" API int map_fit(
//     int map_id, int data_id,
//     int* status
// );