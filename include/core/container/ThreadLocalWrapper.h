// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/GeneralSettings.h>
#include <utility/MultiThreading.h>

#include <vector>
#include <functional>
#include <unordered_map>

namespace ausaxs::container {
    /**
     * @brief A simple wrapper around T to keep track of the thread-local instances of T.
     *        This allows access to all the thread-local data from any single thread.
     *        Note that it is assumed that all threads have a longer lifetime than this class.
     *        
     *        ! Making a static instance of this class will cause Windows DLL to deadlock when the program is closed.
     *        ? This is probably due to the threads owning the data no longer existing when this class is destroyed, combined with the FreeLibrary locking the system resources necessary for C++ to solve this. 
     */
    template <typename T>
    class ThreadLocalWrapper {
        public:
            /**
             * @brief Create a wrapper around T, and create a thread-local instance of T for each thread using the given arguments.
             * 
             * ! This constructor *must* be called from the main thread, otherwise it will not receive its own thread-local instance.
             */
            template <typename... Args>
            ThreadLocalWrapper(Args&&... args) {
                auto ids = utility::multi_threading::get_global_pool()->get_thread_ids();
                for (auto& id : ids) {
                    data.emplace(id, T(args...));
                }
                data.emplace(std::this_thread::get_id(), T(std::move(args)...));
            }

            /**
             * @brief Get the thread-local instance of the wrapped type.
             */
            T& get() {return data.at(std::this_thread::get_id());}

            // @copydoc get()
            const T& get() const {return data.at(std::this_thread::get_id());}

            /**
             * @brief Get the thread-local instances of the wrapped type for all threads.
             */
            std::vector<std::reference_wrapper<T>> get_all() {
                std::vector<std::reference_wrapper<T>> result; result.reserve(data.size());
                for (auto& [id, t] : data) {result.emplace_back(t);}
                return result;
            }

            /**
             * @brief Get the number of thread-local instances of the wrapped type.
             */
            std::size_t size() const {return data.size();}

            /**
             * @brief Reinitialize all thread-local instances of the wrapped type using the given arguments.
             */
            template <typename... Args>
            void reinitialize_all(Args&&... args) {
                for (auto& e : data) {
                    e.second = T(args...);
                }
            }

            /**
             * @brief Merge all thread-local instances of the wrapped type into a single instance.
             */
            T merge() const {
                T result = get();
                auto this_id = std::this_thread::get_id();
                if constexpr (std::ranges::range<T>) {
                    for (const auto& [id, t] : data) {
                        if (id == this_id) {continue;}
                        std::transform(t.begin(), t.end(), result.begin(), result.begin(), std::plus<>());
                    }
                    return result;
                } else {
                    for (const auto& [id, t] : data) {
                        if (id == this_id) {continue;}
                        result += t;
                    }
                    return result;
                }
            }

        private:
            std::unordered_map<std::thread::id, T> data;
    };
}