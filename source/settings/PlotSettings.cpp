#include <settings/PlotSettings.h>
#include <settings/SettingsIORegistry.h>

namespace settings::plots {
    std::string format = "png";
    std::vector<double> contour = {};

    namespace io {
        settings::io::SettingSection plot_settings("Plot", {
            settings::io::create(format, "format"),
            settings::io::create(contour, "contour")
        });
    }
}