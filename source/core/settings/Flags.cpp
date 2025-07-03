// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/Flags.h>

using namespace ausaxs::settings;

bool flags::data_rebin = false;
char flags::last_parsed_unit = ' ';
bool flags::init_histogram_manager = true;