#pragma once

#include <utility/observer_ptr.h>

#include <vector>
#include <memory>
#include <concepts>

namespace utility {
    namespace detail {
        template<typename T> concept observer_type = requires(T t, const T& value) {
            { t.notify(value) } -> std::same_as<void>;
        };
    }

    template<typename T, detail::observer_type Observer>
    class Observable {
        using observer_ptr = std::observer_ptr<Observer>;
        public: 
            Observable() = default;

            void attach_observer(observer_ptr o) {
                observers.push_back(o);
            }

            void detach_observer(observer_ptr o) {
                observers.erase(std::remove(observers.begin(), observers.end(), &o), observers.end());
            }

            void notify(const T& value) {
                for (auto& o : observers) {
                    o->notify(value);
                }
            }

        private:
            std::vector<observer_ptr> observers;
    };
}