#include <utility/MultiThreading.h>
#include <settings/GeneralSettings.h>

observer_ptr<BS::thread_pool> utility::multi_threading::get_global_pool() {
    // statics in functions are initialized on first call, so ok to use settings::general::threads here
    static std::unique_ptr<BS::thread_pool> pool = std::make_unique<BS::thread_pool>(settings::general::threads);
    return pool.get();
}