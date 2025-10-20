// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/SimpleExvModel.h>
#include <utility/Exceptions.h>
#include <constants/Constants.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/ExvSettings.h>
#include <settings/FitSettings.h>
#include <settings/SettingsIORegistry.h>

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

settings::io::SettingSection hist_section("Histogram", {
    settings::io::create(settings::hist::weighted_bins, "weighted_bins")
});

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