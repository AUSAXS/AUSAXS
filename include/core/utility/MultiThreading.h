#pragma once

#include <utility/observer_ptr.h>

#include <BS_thread_pool.hpp>

namespace ausaxs::utility::multi_threading {
    /**
     * @brief Get the global thread pool.
     *        This pool is initialized upon first call, so make sure to set the number of threads (settings::general::threads) before calling this function.
     */
    observer_ptr<BS::light_thread_pool> get_global_pool();
}