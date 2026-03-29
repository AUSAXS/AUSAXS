// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/lookup/FormFactorProduct.h>
#include <container/ArrayContainer2D.h>

namespace ausaxs::form_factor {
    namespace lookup::exv {using table_t = container::ArrayContainer2D<FormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>;}
    namespace lookup::cross {using table_t = container::ArrayContainer2D<FormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>;}
    namespace lookup::atomic {using table_t = container::ArrayContainer2D<FormFactorProduct, form_factor::get_count(), form_factor::get_count()>;}
}