/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/MoleculeSettings.h>
#include <settings/SettingsIORegistry.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

bool settings::molecule::center = true;
bool settings::molecule::implicit_hydrogens = true;
bool settings::molecule::use_occupancy = true;

settings::molecule::DisplacedVolumeSet settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Default;

#if DEBUG
    bool settings::molecule::throw_on_unknown_atom = true;
#else
    bool settings::molecule::throw_on_unknown_atom = false;
#endif

namespace ausaxs::settings::molecule::io {
    settings::io::SettingSection molecule_settings("Molecule", {
        settings::io::create(center, "center"),
        settings::io::create(throw_on_unknown_atom, "throw_on_unknown_atom"),
        settings::io::create(implicit_hydrogens, "implicit_hydrogens"),
        settings::io::create(use_occupancy, "use_occupancy"),
        settings::io::create(displaced_volume_set, "exv_volume")
    });

    settings::io::SettingSection hydrate_settings("Hydrate", {
        settings::io::create(hydrate::hydration_strategy, "hydration_strategy"),
        settings::io::create(hydrate::culling_strategy, "culling_strategy")
    });
}

settings::hydrate::HydrationStrategy settings::hydrate::hydration_strategy = HydrationStrategy::RadialStrategy;
settings::hydrate::CullingStrategy settings::hydrate::culling_strategy = CullingStrategy::NoStrategy;
double settings::hydrate::shell_correction = -0.35;

template<> std::string settings::io::detail::SettingRef<settings::hydrate::HydrationStrategy>::get() const {
    switch (settingref) {
        case settings::hydrate::HydrationStrategy::RadialStrategy: return "radial";
        case settings::hydrate::HydrationStrategy::AxesStrategy: return "axes";
        case settings::hydrate::HydrationStrategy::JanStrategy: return "jan";
        case settings::hydrate::HydrationStrategy::PepsiStrategy: return "pepsi";
        case settings::hydrate::HydrationStrategy::NoStrategy: return "none";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::hydrate::HydrationStrategy>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]); 
    if (str == "radial") {settingref = settings::hydrate::HydrationStrategy::RadialStrategy;}
    else if (str == "axes") {settingref = settings::hydrate::HydrationStrategy::AxesStrategy;}
    else if (str == "jan") {settingref = settings::hydrate::HydrationStrategy::JanStrategy;}
    else if (str == "pepsi") {settingref = settings::hydrate::HydrationStrategy::PepsiStrategy;}
    else if (str == "none") {settingref = settings::hydrate::HydrationStrategy::NoStrategy;}
    else if (!val[0].empty() && std::isdigit(val[0][0])) {settingref = static_cast<settings::hydrate::HydrationStrategy>(std::stoi(val[0]));}
    else {
        throw except::io_error("settings::hydrate::placement_strategy: Unkown HydrationStrategy. Did you forget to add parsing support for it in MoleculeSettings.cpp?");
    }
}

template<> std::string settings::io::detail::SettingRef<settings::hydrate::CullingStrategy>::get() const {return std::to_string(static_cast<int>(settingref));}
template<> void settings::io::detail::SettingRef<settings::hydrate::CullingStrategy>::set(const std::vector<std::string>& val) {
    settingref = static_cast<settings::hydrate::CullingStrategy>(std::stoi(val[0]));
}

template<> std::string settings::io::detail::SettingRef<settings::molecule::DisplacedVolumeSet>::get() const {
    switch (settingref) {
        case settings::molecule::DisplacedVolumeSet::Traube: return "traube";
        case settings::molecule::DisplacedVolumeSet::Voronoi_implicit_H: return "voronoi_implicit_h";
        case settings::molecule::DisplacedVolumeSet::Voronoi_explicit_H: return "voronoi_explicit_h";
        case settings::molecule::DisplacedVolumeSet::MinimumFluctutation_implicit_H: return "minimum_fluctuation_implicit_h";
        case settings::molecule::DisplacedVolumeSet::MinimumFluctutation_explicit_H: return "minimum_fluctuation_explicit_h";
        case settings::molecule::DisplacedVolumeSet::vdw: return "vdw";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::molecule::DisplacedVolumeSet>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]); 
    if (str == "traube") {settingref = settings::molecule::DisplacedVolumeSet::Traube;}
    else if (str == "voronoi_implicit_h" || str == "voronoi") {settingref = settings::molecule::DisplacedVolumeSet::Voronoi_implicit_H;}
    else if (str == "voronoi_explicit_h") {settingref = settings::molecule::DisplacedVolumeSet::Voronoi_explicit_H;}
    else if (str == "minimum_fluctuation_implicit_h" || str == "mf") {settingref = settings::molecule::DisplacedVolumeSet::MinimumFluctutation_implicit_H;}
    else if (str == "minimum_fluctuation_explicit_h") {settingref = settings::molecule::DisplacedVolumeSet::MinimumFluctutation_explicit_H;}
    else if (str == "vdw") {settingref = settings::molecule::DisplacedVolumeSet::vdw;}
    else if (!val[0].empty() && std::isdigit(val[0][0])) {settingref = static_cast<settings::molecule::DisplacedVolumeSet>(std::stoi(val[0]));}
    else {
        throw except::io_error("settings::molecule::displaced_volume_set: Unkown DisplacedVolumeSet. Did you forget to add parsing support for it in MoleculeSettings.cpp?");
    }
}