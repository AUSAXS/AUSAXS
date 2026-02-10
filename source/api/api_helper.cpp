// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_helper.h>

extern "C" API void get_last_error_msg(
    char** buffer,
    int* status
) {
    *status = 1;
    *buffer = const_cast<char*>(ErrorMessage::last_error.c_str());
    *status = 0;
}