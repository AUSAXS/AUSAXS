#pragma once

#include <utility/observer_ptr.h>
#include <utility/Observer.h>

#include <list>
#include <memory>

namespace ausaxs::utility {
    namespace detail {
        template<typename T> concept observer_type = requires(T t, const T& value) {
            { t.notify(value) } -> std::same_as<void>;
        };
    }

    template<typename T>
    class Observable {
        public: 
            Observable() = default;
            ~Observable() {
                for (auto& o : observers) {
                    o->on_delete();
                    o->on_delete = [] () {}; // clear the on_delete callback to avoid dangling pointer to this
                }
            }

            /**
             * @brief Notify all observers of a new value.
             */
            void notify(const T& value) {
                for (auto& o : observers) {
                    o->notify(value);
                }
            }

            /**
             * @brief Create a new observer and attach it to this observable. 
             * 
             * When the returned observer is destroyed, it will detach itself from this observable.
             */
            std::unique_ptr<Observer<T>> make_observer() {
                auto o = std::make_unique<Observer<T>>();
                o->on_delete = [this, ptr=o.get()] () {
                    this->detach_observer(ptr);
                };
                this->attach_observer(o.get());
                return o;
            }

        private:
            void attach_observer(observer_ptr<Observer<T>> observer) {
                observers.push_back(observer);
            }

            void detach_observer(observer_ptr<Observer<T>> observer) {
                observers.remove(observer);
            }

            std::list<observer_ptr<Observer<T>>> observers;
    };
}