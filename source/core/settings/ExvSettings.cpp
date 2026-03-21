// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/ExvSettings.h>
#include <settings/SettingsIORegistry.h>
#include <utility/StringUtils.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

settings::exv::ExvMethod settings::exv::exv_method = settings::exv::ExvMethod::Simple;
settings::exv::ExvSet settings::exv::exv_set = settings::exv::ExvSet::Default;
settings::io::SettingSection exv_section("Excluded volume", {
    settings::io::create(settings::exv::exv_method, "exv_model"),
    settings::io::create(settings::exv::exv_set, "exv_volume")
});

template<> std::string settings::io::detail::SettingRef<settings::exv::ExvMethod>::get() const {
    switch (settingref) {
        case settings::exv::ExvMethod::Simple:      return "simple";
        case settings::exv::ExvMethod::Average:     return "average";
        case settings::exv::ExvMethod::Fraser:      return "fraser";
        case settings::exv::ExvMethod::Grid:        return "grid-base";
        case settings::exv::ExvMethod::GridScalable:return "grid-scalable";
        case settings::exv::ExvMethod::GridSurface: return "grid";
        case settings::exv::ExvMethod::CRYSOL:      return "crysol";
        case settings::exv::ExvMethod::FoXS:        return "foxs";
        case settings::exv::ExvMethod::Pepsi:       return "pepsi";
        case settings::exv::ExvMethod::WAXSiS:      return "waxsis";
        case settings::exv::ExvMethod::None:        return "none";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::exv::ExvMethod>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]);
    if (     str == "simple") {settingref = settings::exv::ExvMethod::Simple;}
    else if (str == "average") {settingref = settings::exv::ExvMethod::Average;}
    else if (str == "fraser") {settingref = settings::exv::ExvMethod::Fraser;}
    else if (str == "grid-base") {settingref = settings::exv::ExvMethod::Grid;}
    else if (str == "grid-scalable") {settingref = settings::exv::ExvMethod::GridScalable;}
    else if (str == "grid") {settingref = settings::exv::ExvMethod::GridSurface;}
    else if (str == "crysol") {settingref = settings::exv::ExvMethod::CRYSOL;}
    else if (str == "foxs") {settingref = settings::exv::ExvMethod::FoXS;}
    else if (str == "pepsi") {settingref = settings::exv::ExvMethod::Pepsi;}
    else if (str == "waxsis") {settingref = settings::exv::ExvMethod::WAXSiS;}
    else if (str == "none") {settingref = settings::exv::ExvMethod::None;}
    else if (!val[0].empty() && std::isdigit(val[0][0])) {settingref = static_cast<settings::exv::ExvMethod>(std::stoi(val[0]));}
    else {
        throw except::io_error("settings: Unknown excluded volume method \"" + str + "\". Did you forget to add parsing support for it in ExvSettings.cpp?");
    }
}

template<> std::string settings::io::detail::SettingRef<settings::exv::ExvSet>::get() const {
    switch (settingref) {
        case settings::exv::ExvSet::Traube: return "traube";
        case settings::exv::ExvSet::Voronoi_implicit_H: return "voronoi_implicit_h";
        case settings::exv::ExvSet::Voronoi_explicit_H: return "voronoi_explicit_h";
        case settings::exv::ExvSet::MinimumFluctutation_implicit_H: return "minimum_fluctuation_implicit_h";
        case settings::exv::ExvSet::MinimumFluctutation_explicit_H: return "minimum_fluctuation_explicit_h";
        case settings::exv::ExvSet::vdw: return "vdw";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::exv::ExvSet>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]); 
    if (str == "traube") {settingref = settings::exv::ExvSet::Traube;}
    else if (str == "voronoi_implicit_h" || str == "voronoi") {settingref = settings::exv::ExvSet::Voronoi_implicit_H;}
    else if (str == "voronoi_explicit_h") {settingref = settings::exv::ExvSet::Voronoi_explicit_H;}
    else if (str == "minimum_fluctuation_implicit_h" || str == "mf") {settingref = settings::exv::ExvSet::MinimumFluctutation_implicit_H;}
    else if (str == "minimum_fluctuation_explicit_h") {settingref = settings::exv::ExvSet::MinimumFluctutation_explicit_H;}
    else if (str == "vdw") {settingref = settings::exv::ExvSet::vdw;}
    else if (!val[0].empty() && std::isdigit(val[0][0])) {settingref = static_cast<settings::exv::ExvSet>(std::stoi(val[0]));}
    else {
        throw except::io_error("settings::molecule::displaced_volume_set: Unkown DisplacedVolumeSet \"" + str + "\". Did you forget to add parsing support for it in MoleculeSettings.cpp?");
    }
}