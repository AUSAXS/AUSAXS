// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API bool io_is_pdb(
    const char* path,
    int* status
);

extern "C" API bool io_is_saxs_data(
    const char* path,
    int* status
);

extern "C" API bool io_is_em_map(
    const char* path,
    int* status
);

extern "C" API bool io_is_rigidbody_config(
    const char* path,
    int* status
);