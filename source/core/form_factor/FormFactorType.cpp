// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactor.h>

using namespace ausaxs::form_factor;

double ausaxs::constants::charge::get_ff_charge(form_factor_t type) {
    return lookup::atomic::raw::get(type).I0();
}