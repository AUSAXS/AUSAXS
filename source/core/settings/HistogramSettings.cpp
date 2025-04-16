/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/HistogramSettings.h>
#include <settings/SettingsIORegistry.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>
#include <constants/Constants.h>

using namespace ausaxs;

double settings::axes::qmin = constants::axes::q_axis.min;
double settings::axes::qmax = 0.5;
unsigned int settings::axes::skip = 0;
bool settings::hist::weighted_bins = true;

namespace ausaxs::settings::io {
    settings::io::SettingSection axes_section("Axes", {
        settings::io::create(axes::skip, "skip"),
        settings::io::create(axes::qmin, "qmin"),
        settings::io::create(axes::qmax, "qmax"),
    });
}

settings::hist::HistogramManagerChoice settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PartialHistogramManagerMT;
settings::io::SettingSection hist_section("Histogram", {
    settings::io::create(settings::hist::histogram_manager, "histogram_manager"),
    settings::io::create(settings::hist::weighted_bins, "weighted_bins")
});

template<> std::string settings::io::detail::SettingRef<settings::hist::HistogramManagerChoice>::get() const {
    switch (settingref) {
        case settings::hist::HistogramManagerChoice::HistogramManager: return "hm";
        case settings::hist::HistogramManagerChoice::HistogramManagerMT: return "hmmt";
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg: return "hmmtff";
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit: return "hmmtffx";
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: return "hmmtffg";
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface: return "hmmtffgs";
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv: return "hmmtffgse";
        case settings::hist::HistogramManagerChoice::PartialHistogramManager: return "phm";
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT: return "phmmt";
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMTFFAvg: return "phmmtff";
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMTFFExplicit: return "phmmtffx";
        case settings::hist::HistogramManagerChoice::FoXSManager: return "foxs";
        case settings::hist::HistogramManagerChoice::PepsiManager: return "pepsi";
        case settings::hist::HistogramManagerChoice::CrysolManager: return "crysol";
        // case settings::hist::HistogramManagerChoice::DebugManager: return "debug";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::hist::HistogramManagerChoice>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]);
    if (     str == "hm") {settingref = settings::hist::HistogramManagerChoice::HistogramManager;}
    else if (str == "hmmt" || str == "simple") {settingref = settings::hist::HistogramManagerChoice::HistogramManagerMT;}
    else if (str == "hmmtff") {settingref = settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg;}
    else if (str == "hmmtffx" || str == "fraser") {settingref = settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;}
    else if (str == "hmmtffg") {settingref = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;}
    else if (str == "hmmtffgs" || str == "grid") {settingref = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface;}
    else if (str == "hmmtffgse" || str == "scalable") {settingref = settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv;}
    else if (str == "phm") {settingref = settings::hist::HistogramManagerChoice::PartialHistogramManager;}
    else if (str == "phmmt") {settingref = settings::hist::HistogramManagerChoice::PartialHistogramManagerMT;}
    else if (str == "phmmtff") {settingref = settings::hist::HistogramManagerChoice::PartialHistogramManagerMTFFAvg;}
    else if (str == "phmmtffx") {settingref = settings::hist::HistogramManagerChoice::PartialHistogramManagerMTFFExplicit;}
    // else if (str == "debug") {settingref = settings::hist::HistogramManagerChoice::DebugManager;}
    else if (str == "foxs") {settingref = settings::hist::HistogramManagerChoice::FoXSManager;}
    else if (str == "pepsi") {settingref = settings::hist::HistogramManagerChoice::PepsiManager;}
    else if (str == "crysol") {settingref = settings::hist::HistogramManagerChoice::CrysolManager;}
    else if (!val[0].empty() && std::isdigit(val[0][0])) {settingref = static_cast<settings::hist::HistogramManagerChoice>(std::stoi(val[0]));}
    else {
        throw except::io_error("settings::hist::histogram_manager: Unkown HistogramManagerChoice. Did you forget to add parsing support for it in HistogramSettings.cpp?");
    }
}