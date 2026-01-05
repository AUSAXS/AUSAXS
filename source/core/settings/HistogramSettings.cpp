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
bool settings::hist::weighted_bins = true;

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

auto small_d_range_warning = [] (double bin_width, unsigned int bin_count) {
    static bool warned = false;
    if (warned) {return;}
    warned = true;
    console::print_warning(
        "settings::axes::bin_width: The specified bin width (" + std::to_string(bin_width) + "\u212B) and bin count (" + std::to_string(bin_count) + ") "
        "result in a maximum d-value of less than the recommended " + std::to_string(int(constants::axes::d_axis.max)) + "\u212B. "
        "Structures larger than the new minimum of " + std::to_string(int(bin_width*bin_count)) + "\u212B will trigger segmentation faults. "
        "Consider increasing the bin count or bin width to cover a longer range."
    );
};

// bin_width
settings::detail::Setting<double> settings::axes::bin_width = {
    constants::axes::d_axis.width(),
    [](double& new_width) {
        if (new_width*settings::axes::bin_count < constants::axes::d_axis.max) {small_d_range_warning(new_width, settings::axes::bin_count);}
        if (std::abs(constants::axes::d_axis.width() - new_width) < 1e-6) {
            settings::flags::custom_bin_width = false;
        } else {
            settings::flags::custom_bin_width = true;
        }
        settings::flags::inv_bin_width = 1./new_width;
    }
};

// bin_count
settings::detail::Setting<unsigned int> settings::axes::bin_count = {
    constants::axes::d_axis.bins,
    [](unsigned int& new_count) {
        if (new_count*settings::axes::bin_width < constants::axes::d_axis.max) {small_d_range_warning(settings::axes::bin_width, new_count);}
    }
};

namespace ausaxs::settings::io {
    settings::io::SettingSection axes_section("Axes", {
        settings::io::create(axes::skip, "skip"),
        settings::io::create(axes::qmin, "qmin"),
        settings::io::create(axes::qmax, "qmax"),
    });

    settings::io::SettingSection hist_section("Histogram", {
        settings::io::create(settings::hist::weighted_bins, "weighted_bins"),
        settings::io::create(settings::axes::bin_width, "bin_width"),
        settings::io::create(settings::axes::bin_count, "bin_count"),
    });
}

settings::hist::HistogramManagerChoice settings::hist::get_histogram_manager() {
    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::Simple:
            // if no multi-threading is enabled, switch to the single-threaded manager
            return settings::general::threads == 1 
                ? settings::hist::HistogramManagerChoice::HistogramManager 
                : settings::hist::HistogramManagerChoice::HistogramManagerMT;

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
            return settings::general::threads == 1 
                ? settings::hist::HistogramManagerChoice::HistogramManager 
                : settings::hist::HistogramManagerChoice::HistogramManagerMT;

        default:
            throw except::unexpected("settings::hist::get_histogram_manager: Unknown ExvMethod. Did you forget to add it to the switch statement?");
    }
}