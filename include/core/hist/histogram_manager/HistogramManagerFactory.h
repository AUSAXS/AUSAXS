// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <utility/observer_ptr.h>
#include <hist/HistFwd.h>
#include <settings/HistogramSettings.h>
#include <settings/Flags.h>

#include <memory>

namespace ausaxs::hist {
    namespace factory {
        std::unique_ptr<IHistogramManager> construct_histogram_manager(
            observer_ptr<const data::Molecule> protein, 
            bool weighted_bins = settings::hist::weighted_bins, 
            bool variable_bin_width = settings::flags::custom_bin_width
        );

        std::unique_ptr<IHistogramManager> construct_histogram_manager(
            observer_ptr<const data::Molecule> protein, 
            settings::hist::HistogramManagerChoice choice, 
            bool weighted_bins = settings::hist::weighted_bins,
            bool variable_bin_width = settings::flags::custom_bin_width
        );
    }
}