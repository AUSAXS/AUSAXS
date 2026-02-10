// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API void io_is_pdb(
    const char* script,
    int* status
);

extern "C" API bool io_is_saxs_data(
    const char* script,
    int* status
);

extern "C" API bool io_is_em_map(
    const char* script,
    int* status
);

extern "C" API bool io_is_rigidbody_config(
    const char* script,
    int* status
);