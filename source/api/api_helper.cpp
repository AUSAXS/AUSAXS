// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_helper.h>

#include <iostream>

void report_api_exception(const char* what) noexcept {
    try {
        std::cerr << "AUSAXS error: " << (what ? what : "An unknown error occurred.") << std::endl;
    } catch (...) {} // reporting a failure must never itself throw out of the handler that called it
}

extern "C" API void get_last_error_msg(
    char** buffer,
    int* status
) {
    *status = 1;
    *buffer = const_cast<char*>(ErrorMessage::last_error.c_str());
    *status = 0;
}