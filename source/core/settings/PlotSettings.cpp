/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/PlotSettings.h>
#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

std::string settings::plots::format = "png";
std::vector<double> settings::plots::contour = {};

namespace ausaxs::settings::io {
    settings::io::SettingSection plot_section("Plot", {
        settings::io::create(plots::format, "format")
    });
}