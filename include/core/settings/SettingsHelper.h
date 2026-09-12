// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <functional>

namespace ausaxs::settings::detail {
    template<typename T>
    struct Setting {
        T value;
        std::function<void(T&)> on_change;
        // Returns the stored value rather than *this, so that the result of `setting = x` can be bound as a T&;
        // NOLINTNEXTLINE(misc-unconventional-assign-operator)
        T& operator=(T other) {
            value = std::move(other);
            if (on_change) {on_change(value);}
            return value;
        }
        operator T() const {return value;}
    };
}