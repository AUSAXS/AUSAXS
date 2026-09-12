// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ArrayContainer2D.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorProduct.h>

namespace ausaxs::form_factor::lookup {
    using table_t = container::ArrayContainer2D<FormFactorProduct, form_factor::total_ff_count, form_factor::total_ff_count>;
}