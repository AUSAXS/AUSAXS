/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/MDSettings.h>
#include <settings//SettingsIORegistry.h>

using namespace ausaxs;

std::string settings::md::gmx_path = "gmx";
std::string settings::md::gmx_source_path = "";
std::string settings::md::buffer_path = "";

std::string settings::md::gmx_top_path() {
    return gmx_source_path + "/share/top/";
}

namespace ausaxs::settings::md::io {
    settings::io::SettingSection md_section("MD", {
        settings::io::create(gmx_path, {"gmx_exe", "gmx_executable", "gmx"}),
        settings::io::create(gmx_source_path, {"gmx_source", "gmx_source_path"}),
        settings::io::create(buffer_path, {"buffer_path", "buffer"}),
    });
}
