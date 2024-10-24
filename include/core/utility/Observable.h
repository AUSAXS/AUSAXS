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
        class Observer {
            friend class Observable;
            public:
                Observer() = default;

                ~Observer() {
                    on_delete();
                }

                std::function<void(const T&)> on_notify;

            private:
                void notify(const T& value) {
                    on_notify(value);
                }

                std::function<void()> on_delete;
        };

        public: 
            Observable() = default;

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
            std::unique_ptr<Observer> make_observer() {
                auto o = std::make_unique<Observer>();
                o->on_delete = [this, &o] () {
                    this->detach_observer(o.get());
                };
                this->attach_observer(o.get());
                return o;
            }

        private:
            void attach_observer(observer_ptr<Observer> observer) {
                observers.push_back(observer);
            }

            void detach_observer(observer_ptr<Observer> observer) {
                observers.remove(observer);
            }

            std::list<observer_ptr<Observer>> observers;
    };
}