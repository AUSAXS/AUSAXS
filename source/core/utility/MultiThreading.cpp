// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/MultiThreading.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;

observer_ptr<BS::light_thread_pool> utility::multi_threading::get_global_pool() {
    // statics in functions are initialized on first call, so ok to use settings::general::threads here
    static std::unique_ptr<BS::light_thread_pool> pool = std::make_unique<BS::light_thread_pool>(settings::general::threads);
    return pool.get();
}