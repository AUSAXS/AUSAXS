// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <container/ArrayContainer2D.h>

namespace ausaxs::form_factor {
    namespace lookup {using table_t = container::ArrayContainer2D<FormFactorProduct, form_factor::total_ff_count, form_factor::total_ff_count>;}
}