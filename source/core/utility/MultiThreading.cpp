/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/MultiThreading.h>
#include <settings/GeneralSettings.h>

observer_ptr<BS::thread_pool> utility::multi_threading::get_global_pool() {
    // statics in functions are initialized on first call, so ok to use settings::general::threads here
    static std::unique_ptr<BS::thread_pool> pool = std::make_unique<BS::thread_pool>(settings::general::threads);
    return pool.get();
}