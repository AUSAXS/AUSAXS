// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/SimpleExvModel.h>
#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <constants/Constants.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/ExvSettings.h>
#include <settings/FitSettings.h>
#include <settings/Flags.h>
#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

unsigned int settings::axes::skip = 0;

// qmin
settings::detail::Setting<double> settings::axes::qmin = {
    constants::axes::q_axis.min,
    [](double& new_qmin) {
        if (new_qmin < 0. || new_qmin > constants::axes::q_axis.max) {
            console::print_warning(
                "settings::axes::qmin: qmin must be in the range "
                "[" + std::to_string(constants::axes::q_axis.min) + ", " + std::to_string(settings::axes::qmax) + "]. "
                "Clamping to closest value."
            );
            new_qmin = std::clamp(new_qmin, constants::axes::q_axis.min, settings::axes::qmax.value);
        }
    }
};

// qmax
settings::detail::Setting<double> settings::axes::qmax = {
    0.5,
    [](double& new_qmax) {
        if (new_qmax < 0. || new_qmax > constants::axes::q_axis.max) {
            console::print_warning(
                "settings::axes::qmax: qmax must be in the range" 
                "[" + std::to_string(settings::axes::qmin) + ", " + std::to_string(constants::axes::q_axis.max) + "]. "
                "Clamping to closest value."
            );
            new_qmax = std::clamp(new_qmax, settings::axes::qmin.value, constants::axes::q_axis.max);
        }
    }
};

// bin_width
settings::detail::Setting<double> settings::axes::bin_width = {
    constants::axes::d_axis.width(),
    [](double& new_width) {
        if (std::abs(constants::axes::d_axis.width() - new_width) < 1e-6) {
            settings::flags::custom_bin_width = false;
        } else {
            settings::flags::custom_bin_width = true;
        }
        settings::flags::inv_bin_width = 1./new_width;
    }
};

bool settings::axes::clamp_to_qrange = true;

namespace ausaxs::settings::io {
    settings::io::SettingSection axes_section("Axes", {
        settings::io::create(axes::skip, "skip"),
        settings::io::create(axes::qmin, "qmin"),
        settings::io::create(axes::qmax, "qmax"),
        settings::io::create(axes::clamp_to_qrange, "clamp_to_q"),
    });

    settings::io::SettingSection hist_section("Histogram", {
        settings::io::create(settings::hist::weighted_bins, "weighted_bins"),
        settings::io::create(settings::axes::bin_width, "bin_width"),
    });
}

namespace {
    settings::hist::HistogramManagerChoice plain_manager() {
        using Choice = settings::hist::HistogramManagerChoice;
        bool st = settings::general::threads == 1; // if no multi-threading is enabled, switch to the single-threaded manager
        if (settings::flags::prefer_partial_manager) {
            return st ? Choice::PartialHistogramManager : Choice::PartialHistogramManagerMT;
        }
        return st ? Choice::HistogramManager : Choice::HistogramManagerMT;
    }
}

settings::hist::HistogramManagerChoice settings::hist::get_histogram_manager() {
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::Simple:
            return plain_manager();

        case settings::exv::ExvMethod::Average: 
            return settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg;

        case settings::exv::ExvMethod::Fraser:
            return settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;

        case settings::exv::ExvMethod::Grid:
        case settings::exv::ExvMethod::WAXSiS:
            return settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;

        case settings::exv::ExvMethod::GridScalable:
            // if no exv fitting is performed, switch to the faster grid manager 
            return settings::fit::fit_excluded_volume 
                ? settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv 
                : settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;

        case settings::exv::ExvMethod::GridSurface:
            return settings::fit::fit_excluded_volume 
                ? settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface 
                : settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;

        case settings::exv::ExvMethod::CRYSOL:
            return settings::hist::HistogramManagerChoice::CrysolManager;

        case settings::exv::ExvMethod::FoXS:
            return settings::hist::HistogramManagerChoice::FoXSManager;
    
        case settings::exv::ExvMethod::Pepsi:
            return settings::hist::HistogramManagerChoice::PepsiManager;

        case settings::exv::ExvMethod::None:
            ausaxs::hist::detail::SimpleExvModel::disable();
            return plain_manager();

        default:
            throw except::unexpected("settings::hist::get_histogram_manager: Unknown ExvMethod. Did you forget to add it to the switch statement?");
    }
}

bool settings::hist::supports_partial_calculation(settings::hist::HistogramManagerChoice choice) {
    switch (choice) {
        case settings::hist::HistogramManagerChoice::PartialHistogramManager:
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
        case settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT:
            return true;
        default:
            return false;
    }
}

settings::hist::WeightedBins settings::hist::weighted_bins(settings::hist::WeightedBins::Value::True);
settings::hist::WeightedBins::WeightedBins(std::string_view str) {
    if (auto lc = utility::to_lowercase(str); lc == "auto") {value = WeightedBins::Value::Auto;}
    else {value = utility::parse_bool(str) ? WeightedBins::Value::True : WeightedBins::Value::False;}
}