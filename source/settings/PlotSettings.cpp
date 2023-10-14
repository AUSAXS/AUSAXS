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