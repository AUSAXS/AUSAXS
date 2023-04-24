#include <settings/GridSettings.h>
#include <utility/Exceptions.h>
#include <settings/SettingsIORegistry.h>

namespace settings::grid {
    double percent_water = 0.1;
    double ra = 2.4;
    double rh = 1.5;
    double ra_effective = 2.4;
    double width = 1;
    double scaling = 0.25;
    bool cubic = false;
    Limit3D axes = Limit3D(-250, 250, -250, 250, -250, 250);
    PlacementStrategy placement_strategy = PlacementStrategy::RadialStrategy;
    CullingStrategy culling_strategy = CullingStrategy::CounterStrategy;

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

template<> std::string settings::io::detail::SettingRef<Limit3D>::get() const {
    return std::to_string(settingref.x.min) + " " + std::to_string(settingref.x.max) + " "
         + std::to_string(settingref.y.min) + " " + std::to_string(settingref.y.max) + " " 
         + std::to_string(settingref.z.min) + " " + std::to_string(settingref.z.max);
}
template<> void settings::io::detail::SettingRef<Limit3D>::set(const std::vector<std::string>& val) {
    if (val.size() != 6) throw except::io_error("settings::grid::axes: Expected 6 values, got " + std::to_string(val.size()) + ".");
    settingref = Limit3D(std::stod(val[0]), std::stod(val[1]), std::stod(val[2]), std::stod(val[3]), std::stod(val[4]), std::stod(val[5]));
}

template<> std::string settings::io::detail::SettingRef<settings::grid::PlacementStrategy>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::grid::PlacementStrategy>::set(const std::vector<std::string>& val) {
    if (val[0] == "Radial") {settingref = settings::grid::PlacementStrategy::RadialStrategy;}
    else if (val[0] == "Axes") {settingref = settings::grid::PlacementStrategy::AxesStrategy;}
    else if (val[0] == "Jan") {settingref = settings::grid::PlacementStrategy::JanStrategy;}
    else if (!val[0].empty() && std::isdigit(val[0][0])) {settingref = static_cast<settings::grid::PlacementStrategy>(std::stoi(val[3]));}
    else {
        throw except::io_error("settings::grid::placement_strategy: Unkown PlacementStrategy. Did you forget to add parsing support for it in GridSettings.cpp?");
    }
}

template<> std::string settings::io::detail::SettingRef<settings::grid::CullingStrategy>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::grid::CullingStrategy>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::grid::CullingStrategy>(std::stoi(val[0]));
}