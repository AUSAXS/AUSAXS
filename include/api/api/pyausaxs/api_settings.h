// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API int get_setting(
    const char* name,
    const char** value,
    const char** type,
    int* status
);

extern "C" API void set_setting(
    const char* name,
    const char* value,
    int* status
);