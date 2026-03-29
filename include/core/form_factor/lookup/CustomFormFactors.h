// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>

#include <vector>

namespace ausaxs::form_factor::lookup {
    /**
     * @brief Declare a custom set of form factors to be used in the calculations. 
     * @param custom_form_factors The set of form factors and their indices. 
     */
    void set_custom_form_factors(std::vector<form_factor::form_factor_t> custom_form_factors);
}