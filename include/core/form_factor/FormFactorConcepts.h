// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>

#include <concepts>

namespace ausaxs {
    /**
     * @brief A form factor which can be evaluated at a given q value.
     */
    template<typename T>
    concept FormFactorType = requires(const T& t, double q) {
        {t.evaluate(q)} -> std::convertible_to<double>;
    };

    /**
     * @brief A lookup providing a form factor for each form factor type through a static `get(form_factor_t)` method.
     */
    template<typename T>
    concept FormFactorLookupType = requires(form_factor::form_factor_t type) {
        {T::get(type)} -> FormFactorType;
    };
}
