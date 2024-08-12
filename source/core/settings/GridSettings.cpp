/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/GridSettings.h>
#include <settings/SettingsIORegistry.h>
#include <utility/Exceptions.h>

double settings::grid::water_scaling = 0.01;
double settings::grid::cell_width = 1;
double settings::grid::scaling = 0.25;
bool settings::grid::cubic = false;
double settings::grid::min_exv_radius = 2.15;
unsigned int settings::grid::min_bins = 0;

double settings::grid::exv::surface_thickness = 1;
double settings::grid::exv::radius = 0.5;
bool settings::grid::exv::save = false;

namespace settings::grid::detail {
    double min_score = 0.1;
}

namespace settings::grid::io {
    settings::io::SettingSection grid_settings("Grid", {
        settings::io::create(water_scaling, "water_scaling"),
        settings::io::create(cell_width, "width"),
        settings::io::create(scaling, "scaling"),
        settings::io::create(cubic, "cubic"),
        settings::io::create(min_exv_radius, "rvol"),
        settings::io::create(exv::radius, "exv_radius"),
        settings::io::create(exv::save, "save_exv"),
        settings::io::create(exv::surface_thickness, "surface_thickness"),
        settings::io::create(detail::min_score, "detail.min_score"),
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