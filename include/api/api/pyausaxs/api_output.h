// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API typedef void (*ausaxs_output_cb)(const char* str, int len);

extern "C" API void set_output_callback(ausaxs_output_cb cb);

extern "C" API void reset_output_callback();