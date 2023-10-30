#pragma once

#include <utility/view_ptr.h>

#include <BS_thread_pool.hpp>

namespace utility::multi_threading {
    /**
     * @brief Get the global thread pool.
     *        This pool is initialized upon first call, so make sure to set the number of threads (settings::general::threads) before calling this function.
     */
    view_ptr<BS::thread_pool> get_global_pool();
}