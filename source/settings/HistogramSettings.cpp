#include <settings/HistogramSettings.h>
#include <settings/SettingsIORegistry.h>

namespace settings::axes {
    unsigned int max_distance = 2000;
    double distance_bin_width = 1;
    unsigned int bins = 1000;
    double qmin = 1e-4;
    double qmax = 0.5;
    unsigned int skip = 0;

    settings::io::SettingSection axes_settings("Axes", {
        settings::io::create(max_distance, "max_distance"),
        settings::io::create(distance_bin_width, "distance_bin_width"),
        settings::io::create(bins, "bins"),
        settings::io::create(qmin, "qmin"),
        settings::io::create(qmax, "qmax"),
        settings::io::create(skip, "skip")
    });
}