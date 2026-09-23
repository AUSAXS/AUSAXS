// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <functional>
#include <utility>

namespace ausaxs::settings::detail {
    template<typename T>
    struct Setting {
        explicit Setting(T value, std::function<void(T&)> on_change = nullptr) : value(std::move(value)), on_change(std::move(on_change)) {}

        // Settings are global state. A copy would carry on_change along, and assigning it back would silently skip the callback.
        // Save and restore the underlying value instead: `auto old = setting.value; ...; setting = old;`
        Setting(const Setting&) = delete;
        Setting& operator=(const Setting&) = delete;

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
