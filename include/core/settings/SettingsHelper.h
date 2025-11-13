// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <functional>

namespace ausaxs::settings::detail {
    template<typename T>
    struct Setting {
        T value;
        std::function<void(T&)> on_change;
        T& operator=(const T& other) {
            value = other;
            if (on_change) {on_change(value);}
            return value;
        }
        operator T() const {return value;}
    };
}