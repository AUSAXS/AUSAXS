#include <hist/HistogramSettings.h>

namespace settings::axes {
    settings::detail::SmartOption<unsigned int> max_distance(2000, "max-distance");
    settings::detail::SmartOption<double> distance_bin_width(1, "distance-bin-width");
    settings::detail::SmartOption<unsigned int> bins(1000, "bins");
    settings::detail::SmartOption<double> qmin(1e-4, "qmin");
    settings::detail::SmartOption<double> qmax(0.5, "qmax");
    settings::detail::SmartOption<unsigned int> skip(0, "skip");
}