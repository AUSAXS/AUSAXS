// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <functional>
#include <utility>

namespace ausaxs::settings::detail {
    template<typename T>
    struct Setting {
        explicit Setting(T value, std::function<void(T&)> on_change = nullptr) 
            : value(value), default_value(std::move(value)), on_change(std::move(on_change)) 
        {}

        // Settings are global state. A copy would carry on_change along, and assigning it back would silently skip the callback.
        Setting(const Setting&) = delete;
        Setting& operator=(const Setting&) = delete;

        T value;
        const T default_value;
        std::function<void(T&)> on_change;
        // Returns the stored value rather than *this, so that the result of `setting = x` can be bound as a T&;
        // NOLINTNEXTLINE(misc-unconventional-assign-operator)
        T& operator=(T other) {
            value = std::move(other);
            if (on_change) {on_change(value);}
            return value;
        }
        operator const T&() const {return value;}
        const T* operator->() const {return &value;}

        void reset() {*this = default_value;}
    };

    template<typename T, typename U> auto operator+(const Setting<T>& lhs, const U& rhs) -> decltype(lhs.value + rhs) {return lhs.value + rhs;}
    template<typename T, typename U> auto operator+(const U& lhs, const Setting<T>& rhs) -> decltype(lhs + rhs.value) {return lhs + rhs.value;}
    template<typename T, typename U> auto operator+=(Setting<T>& lhs, const U& rhs) -> decltype(lhs.value += rhs) {lhs = lhs.value + rhs; return lhs.value;}
    template<typename T, typename U> auto operator==(const Setting<T>& lhs, const U& rhs) -> decltype(lhs.value == rhs) {return lhs.value == rhs;}
}

namespace ausaxs::settings {
    template<typename T> using Setting = detail::Setting<T>;
}
