// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/MultiThreading.h>

#include <settings/GeneralSettings.h>

using namespace ausaxs;

observer_ptr<utility::multi_threading::ThreadPool> utility::multi_threading::get_global_pool() {
#ifdef AUSAXS_SINGLE_THREADED
    static ThreadPool pool;
    return &pool;
#else
    // statics in functions are initialized on first call, so ok to use settings::general::threads here
    static std::unique_ptr<ThreadPool> pool = std::make_unique<ThreadPool>(settings::general::threads);
    return pool.get();
#endif
}