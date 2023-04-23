#include <hydrate/GridSettings.h>
#include <utility/Exceptions.h>
#include <utility/settings/SettingsIORegistry.h>

namespace settings::grid {
    double percent_water = 0.1;
    double ra = 2.4;
    double rh = 1.5;
    double ra_effective = 2.4;
    double width = 1;
    double scaling = 0.25;
    bool cubic = false;

    namespace detail {
        double min_score = 0.1;
    }

    namespace io {
        settings::io::SettingSection grid_settings("Grid", {
            settings::io::create(percent_water, "percent_water"),
            settings::io::create(ra, "ra"),
            settings::io::create(rh, "rh"),
            settings::io::create(ra_effective, "ra_effective"),
            settings::io::create(width, "width"),
            settings::io::create(scaling, "scaling"),
            settings::io::create(cubic, "cubic"),
            settings::io::create(axes, "axes"),
            settings::io::create(detail::min_score, "detail.min_score")
        });
    }
}

Limit3D settings::grid::axes = Limit3D(-250, 250, -250, 250, -250, 250);
// template<> std::string SmartOption<Limit3D>::get() const {
//     return std::to_string(value.x.min) + " " + std::to_string(value.x.max) + " " + std::to_string(value.y.min) + " " + std::to_string(value.y.max) + " " + std::to_string(value.z.min) + " " + std::to_string(value.z.max);
// }
// template<> void SmartOption<Limit3D>::set(const std::vector<std::string>& val) {
//     if (val.size() != 6) throw except::io_error("settings::grid::axes: Expected 6 values, got " + std::to_string(val.size()) + ".");
//     value = Limit3D(std::stod(val[0]), std::stod(val[1]), std::stod(val[2]), std::stod(val[3]), std::stod(val[4]), std::stod(val[5]));
// }

settings::grid::PlacementStrategy settings::grid::placement_strategy = PlacementStrategy::RadialStrategy;
// template<> std::string SmartOption<settings::grid::PlacementStrategy>::get() const {return std::to_string(static_cast<int>(value));}
// template<> void SmartOption<settings::grid::PlacementStrategy>::set(const std::vector<std::string>& val) {
//     if (val[0] == "Radial") {value = settings::grid::PlacementStrategy::RadialStrategy;}
//     else if (val[0] == "Axes") {value = settings::grid::PlacementStrategy::AxesStrategy;}
//     else if (val[0] == "Jan") {value = settings::grid::PlacementStrategy::JanStrategy;}
//     else if (!val[0].empty() && std::isdigit(val[0][0])) {value = static_cast<settings::grid::PlacementStrategy>(std::stoi(val[3]));}
//     else {
//         throw except::io_error("settings::grid::placement_strategy: Unkown PlacementStrategy. Did you forget to add parsing support for it in GridSettings.cpp?");
//     }
// }

settings::grid::CullingStrategy settings::grid::culling_strategy = settings::grid::CullingStrategy::CounterStrategy;
// template<> std::string SmartOption<settings::grid::CullingStrategy>::get() const {return std::to_string(static_cast<int>(value));}
// template<> void SmartOption<settings::grid::CullingStrategy>::set(const std::vector<std::string>& val) {
//     value = static_cast<settings::grid::CullingStrategy>(std::stoi(val[0]));
// }