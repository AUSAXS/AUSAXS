// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <settings/HistogramSettings.h>
#include <settings/InternalState.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::hist::factory {
    std::unique_ptr<IHistogramManager> construct_histogram_manager(
        observer_ptr<const data::Molecule> protein, 
        bool weighted_bins = settings::hist::weighted_bins.is_true()
    );

    std::unique_ptr<IHistogramManager> construct_histogram_manager(
        observer_ptr<const data::Molecule> protein, 
        settings::hist::HistogramManagerChoice choice, 
        bool weighted_bins = settings::hist::weighted_bins.is_true()
    );
}