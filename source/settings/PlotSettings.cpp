/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/PlotSettings.h>
#include <settings/SettingsIORegistry.h>

std::string settings::plots::format = "png";
std::vector<double> settings::plots::contour = {};

namespace settings::plots::io {
    settings::io::SettingSection plot_settings("Plot", {
        settings::io::create(format, "format"),
        settings::io::create(contour, "contour")
    });
}