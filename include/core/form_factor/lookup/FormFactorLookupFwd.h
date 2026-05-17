// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorProduct.h>
#include <container/ArrayContainer2D.h>
#include <settings/FormFactorSettings.h>

namespace ausaxs::form_factor {
    namespace lookup {using table_t = container::ArrayContainer2D<FormFactorProduct, settings::form_factor::max_ff_types, settings::form_factor::max_ff_types>;}
}