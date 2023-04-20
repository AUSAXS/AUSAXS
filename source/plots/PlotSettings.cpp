#include <plots/PlotSettings.h>

namespace settings::plots {
    settings::detail::SmartOption<std::string> format("png", "format");
    settings::detail::SmartOption<std::vector<double>> contour({}, {"contours", "contour-levels"});
}