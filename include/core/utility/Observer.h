#pragma once

#include <functional>

namespace ausaxs::utility {
    namespace detail {
        template<typename T> concept observable_type = requires(T t, const T& value) {
            { t.attach_observer(value) } -> std::same_as<void>;
            { t.detach_observer(value) } -> std::same_as<void>;
            { t.notify(value) } -> std::same_as<void>;
        };
    }

    template<typename T>
    struct Observer {
        Observer() = default;

        void notify(const T& value) {
            on_notify(value);
        }

        std::function<void(const T&)> on_notify;
    };
}