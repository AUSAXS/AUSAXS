// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/observer_ptr.h>

#ifdef AUSAXS_SINGLE_THREADED
    #include <thread>
    #include <utility>
    #include <vector>
#else
    #include <BS_thread_pool.hpp>
#endif

namespace ausaxs::utility::multi_threading {
#ifdef AUSAXS_SINGLE_THREADED
    class SerialThreadPool {
        public:
            template<typename Task>
            void detach_task(Task&& task) {std::forward<Task>(task)();}

            template<typename Index, typename Task>
            void detach_blocks(Index first, Index last, Task&& task) {
                if (first < last) {std::forward<Task>(task)(first, last);}
            }

            template<typename Index, typename Task>
            void detach_loop(Index first, Index last, Task&& task) {
                for (Index i = first; i < last; ++i) {task(i);}
            }

            struct CompletedTasks {
                void wait() const noexcept {}
            };

            template<typename Index, typename Task>
            CompletedTasks submit_blocks(Index first, Index last, Task&& task) {
                detach_blocks(first, last, std::forward<Task>(task));
                return {};
            }

            void wait() const noexcept {}
            std::vector<std::thread::id> get_thread_ids() const {return {};}
    };
    using ThreadPool = SerialThreadPool;
#else
    using ThreadPool = BS::light_thread_pool;
#endif

    /**
     * @brief Get the global thread pool.
     *        This pool is initialized upon first call, so make sure to set the number of threads (settings::general::threads) before calling this function.
     */
    observer_ptr<ThreadPool> get_global_pool();
}