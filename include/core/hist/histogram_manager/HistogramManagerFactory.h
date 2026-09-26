// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <settings/ExvSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/InternalState.h>
#include <utility/observer_ptr.h>

#include <memory>

namespace ausaxs::hist::factory {
    /**
     * @brief Construct the manager of the kind chosen by the settings, upgraded to its symmetry-aware kind if @a protein has symmetries.
     *        See the other overload for how the variant follows from @a exv_method.
     */
    std::unique_ptr<IHistogramManager> construct_histogram_manager(
        observer_ptr<const data::Molecule> protein,
        bool weighted_bins = settings::hist::weighted_bins.is_true(),
        settings::exv::ExvMethod exv_method = settings::exv::exv_method
    );

    /**
     * @brief Construct a manager of the kind @a choice, in the variant of the excluded volume model @a exv_method.
     *
     * The simple models get the weighted manager of the kind, and the form factor models its form factor-resolved variant.
     * The single-threaded kinds are weighted reference implementations, so the form factor models use their multithreaded kinds.
     * The grid models always get their own manager, which recalculates everything.
     */
    std::unique_ptr<IHistogramManager> construct_histogram_manager(
        observer_ptr<const data::Molecule> protein,
        settings::hist::HistogramManagerChoice choice,
        bool weighted_bins = settings::hist::weighted_bins.is_true(),
        settings::exv::ExvMethod exv_method = settings::exv::exv_method
    );

    /**
     * @brief Whether @a exv_method is one of the grid models, which evaluate the excluded volume on a grid of their own.
     */
    bool uses_grid_exv(settings::exv::ExvMethod exv_method = settings::exv::exv_method);
}