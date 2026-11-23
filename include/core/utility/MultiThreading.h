// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/observer_ptr.h>

#include <functional>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

namespace ausaxs::utility::multi_threading {
    class ThreadPoolImpl;

    struct Tasks {
        void wait() const;
        std::function<void()> wait_;
    };

    class ThreadPool {
        public:
            explicit ThreadPool(int threads);
            ~ThreadPool();

            ThreadPool(const ThreadPool&) = delete;
            ThreadPool& operator=(const ThreadPool&) = delete;

            void detach_task(std::function<void()> task);
            void detach_blocks(int first, int last, std::function<void(int, int)> task);

            /**
             * @brief Run task(i) for every i in [first, last), split into blocks.
             */
            template<typename F>
            void detach_loop(int first, int last, F&& task) {
                detach_blocks(first, last, [task = std::forward<F>(task)] (int start, int end) {
                    for (int i = start; i < end; ++i) {task(i);}
                });
            }

            Tasks submit_blocks(int first, int last, std::function<void(int, int)> task);
            void wait();
            std::vector<std::thread::id> get_thread_ids() const;

        private:
            std::unique_ptr<ThreadPoolImpl> impl;
    };

    /**
     * @brief Get the global thread pool.
     *        This pool is initialized upon first call, so make sure to set the number of threads (settings::general::threads) before calling this function.
     */
    observer_ptr<ThreadPool> get_global_pool();
}
