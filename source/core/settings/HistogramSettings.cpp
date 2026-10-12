// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/HistogramSettings.h>

#include <constants/ConstantsAxes.h>
#include <settings/ExvSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/InternalState.h>
#include <settings/SettingsIORegistry.h>
#include <utility/Console.h>
#include <utility/StringUtils.h>

#include <algorithm>

using namespace ausaxs;

settings::Setting<int> settings::axes::skip{0};

// qmin
settings::Setting<double> settings::axes::qmin{
    constants::axes::q_axis.min,
    [](double& new_qmin) {
        if (new_qmin < 0. || new_qmin > constants::axes::q_axis.max) {
            console::print_warning(
                "settings::axes::qmin: qmin must be in the range "
                "[" + std::to_string(constants::axes::q_axis.min) + ", " + std::to_string(settings::axes::qmax) + "]. "
                "Clamping to closest value."
            );
            new_qmin = std::clamp<double>(new_qmin, constants::axes::q_axis.min, settings::axes::qmax);
        }
    }
};

// qmax
settings::Setting<double> settings::axes::qmax{
    0.5,
    [](double& new_qmax) {
        if (new_qmax < 0. || new_qmax > constants::axes::q_axis.max) {
            console::print_warning(
                "settings::axes::qmax: qmax must be in the range" 
                "[" + std::to_string(settings::axes::qmin) + ", " + std::to_string(constants::axes::q_axis.max) + "]. "
                "Clamping to closest value."
            );
            new_qmax = std::clamp<double>(new_qmax, settings::axes::qmin, constants::axes::q_axis.max);
        }
    }
};

// bin_width
settings::Setting<double> settings::axes::bin_width{
    constants::axes::d_axis.width(),
    [](double& new_width) {
        settings::internal_state::inv_bin_width = 1./new_width;
    }
};

settings::Setting<bool> settings::axes::clamp_to_qrange{true};
settings::Setting<bool> settings::axes::rebin{false};

namespace {
    using namespace ausaxs::settings;
    settings::io::SettingSection axes_section("Axes", {
        settings::io::create(axes::skip, "skip"),
        settings::io::create(axes::qmin, "qmin"),
        settings::io::create(axes::qmax, "qmax"),
        settings::io::create(axes::clamp_to_qrange, "clamp_to_q"),
        settings::io::create(axes::rebin, "rebin"),
    });

    settings::io::SettingSection hist_section("Histogram", {
        settings::io::create(settings::hist::weighted_bins, "weighted_bins"),
        settings::io::create(settings::axes::bin_width, "bin_width"),
    });
}

settings::hist::HistogramManagerChoice settings::hist::get_histogram_manager() {
    using Choice = settings::hist::HistogramManagerChoice;
    bool st = settings::general::threads == 1; // if no multi-threading is enabled, switch to the single-threaded manager
    if (settings::internal_state::prefer_partial_manager) {
        return st ? Choice::PartialHistogramManager : Choice::PartialHistogramManagerMT;
    }
    return st ? Choice::HistogramManager : Choice::HistogramManagerMT;
}

bool settings::hist::supports_partial_calculation(settings::hist::HistogramManagerChoice choice) {
    switch (settings::exv::exv_method) {
        // the grid models are only implemented for the managers recalculating everything
        case settings::exv::ExvMethod::Grid:
        case settings::exv::ExvMethod::GridScalable:
        case settings::exv::ExvMethod::GridSurface:
        case settings::exv::ExvMethod::WAXSiS:
            return false;
        default:
            break;
    }

    switch (choice) {
        case settings::hist::HistogramManagerChoice::PartialHistogramManager:
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
        case settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT:
            return true;
        default:
            return false;
    }
}

settings::Setting<settings::hist::WeightedBins> settings::hist::weighted_bins{settings::hist::WeightedBins::Value::True};
settings::hist::WeightedBins::WeightedBins(std::string_view str) {
    if (auto lc = utility::to_lowercase(str); lc == "auto") {value = WeightedBins::Value::Auto;}
    else {value = utility::parse_bool(str) ? WeightedBins::Value::True : WeightedBins::Value::False;}
}