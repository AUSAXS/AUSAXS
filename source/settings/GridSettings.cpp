#include <settings/GridSettings.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <settings/SettingsIORegistry.h>

double settings::grid::water_scaling = 0.01;
double settings::grid::width = 1;
double settings::grid::scaling = 0.25;
bool settings::grid::cubic = false;
double settings::grid::rvol = 2.15;
int settings::grid::exv_radius = 2;
bool settings::grid::save_exv = false;
settings::grid::PlacementStrategy settings::grid::placement_strategy = PlacementStrategy::RadialStrategy;
settings::grid::CullingStrategy settings::grid::culling_strategy = CullingStrategy::CounterStrategy;

namespace settings::grid::detail {
    double min_score = 0.1;
}

namespace settings::grid::io {
    settings::io::SettingSection grid_settings("Grid", {
        settings::io::create(water_scaling, "water_scaling"),
        settings::io::create(width, "width"),
        settings::io::create(scaling, "scaling"),
        settings::io::create(cubic, "cubic"),
        settings::io::create(rvol, "rvol"),
        settings::io::create(exv_radius, "exv_radius"),
        settings::io::create(save_exv, "save_exv"),
        settings::io::create(detail::min_score, "detail.min_score"),
        settings::io::create(placement_strategy, "placement_strategy"),
        settings::io::create(culling_strategy, "culling_strategy")
    });
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

template<> std::string settings::io::detail::SettingRef<settings::grid::PlacementStrategy>::get() const {
    switch (settingref) {
        case settings::grid::PlacementStrategy::RadialStrategy: return "radial";
        case settings::grid::PlacementStrategy::AxesStrategy: return "axes";
        case settings::grid::PlacementStrategy::JanStrategy: return "jan";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::grid::PlacementStrategy>::set(const std::vector<std::string>& val) {
    if (utility::to_lowercase(val[0]) == "radial") {settingref = settings::grid::PlacementStrategy::RadialStrategy;}
    else if (utility::to_lowercase(val[0]) == "axes") {settingref = settings::grid::PlacementStrategy::AxesStrategy;}
    else if (utility::to_lowercase(val[0]) == "jan") {settingref = settings::grid::PlacementStrategy::JanStrategy;}
    else if (!val[0].empty() && std::isdigit(val[0][0])) {settingref = static_cast<settings::grid::PlacementStrategy>(std::stoi(val[0]));}
    else {
        throw except::io_error("settings::grid::placement_strategy: Unkown PlacementStrategy. Did you forget to add parsing support for it in GridSettings.cpp?");
    }
}

template<> std::string settings::io::detail::SettingRef<settings::grid::CullingStrategy>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::grid::CullingStrategy>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::grid::CullingStrategy>(std::stoi(val[0]));
}