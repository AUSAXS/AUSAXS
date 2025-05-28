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
        ~Observer() {
            on_delete();
        }

        std::function<void(const T&)> on_notify;
        std::function<void()> on_delete;

        void notify(const T& value) {
            on_notify(value);
        }
    };

}