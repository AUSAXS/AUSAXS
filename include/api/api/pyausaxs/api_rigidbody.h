// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API int rigidbody_load_script(
    const char* script,
    int* status
);

extern "C" API int rigidbody_validate(
    int rigidbody_id,
    int* status
);

extern "C" API int rigidbody_run(
    int rigidbody_id,
    int* status
);