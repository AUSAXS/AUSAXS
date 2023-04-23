#include <plots/PlotSettings.h>
#include <utility/settings/SettingsIORegistry.h>

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