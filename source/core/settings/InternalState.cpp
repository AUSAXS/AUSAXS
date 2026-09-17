// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/InternalState.h>

#include <constants/ConstantsAxes.h>

using namespace ausaxs::settings;

bool internal_state::custom_bin_width = false;
double internal_state::inv_bin_width = 1./constants::axes::d_axis.width();
bool internal_state::prefer_partial_manager = false;
bool internal_state::allow_decorrelate_atom_order = true;
