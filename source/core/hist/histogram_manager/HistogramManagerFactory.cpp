// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerFactory.h>

#include <data/Molecule.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <settings/HistogramSettings.h>
#include <settings/InternalState.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::hist::factory;

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, bool weighted_bins
) {
    auto choice = settings::hist::get_histogram_manager();
    if (settings::internal_state::prefer_partial_manager && !settings::hist::supports_partial_calculation(choice)) {
        console::print_warning(
            "construct_histogram_manager: A partial histogram manager was requested, but the chosen excluded volume method has no partial implementation. "
            "Every update will recalculate the full histogram."
        );
    }

    if (protein->symmetry().has_symmetries()) {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager:
            case settings::hist::HistogramManagerChoice::HistogramManagerMT:
                choice = settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT;
                break;
            case settings::hist::HistogramManagerChoice::PartialHistogramManager:
            case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
                choice = settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT;
                break;
            default:
                console::print_warning(
                    "construct_histogram_manager: Molecule contains symmetries, but the chosen excluded volume method does not support them. "
                    "Symmetries will be ignored. "
                );
                break;
        }
    }
    return construct_histogram_manager(protein, choice, weighted_bins);
}

namespace {
    template<template<bool> class MANAGER>
    std::unique_ptr<hist::IHistogramManager> create_manager(observer_ptr<const data::Molecule> protein, bool weighted_bins) {
        if (weighted_bins) {
            return std::make_unique<MANAGER<true>>(protein);
        }
        return std::make_unique<MANAGER<false>>(protein);
    }
}

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, settings::hist::HistogramManagerChoice choice, bool weighted_bins
) {
    switch (choice) {
        case settings::hist::HistogramManagerChoice::HistogramManager:
            return create_manager<HistogramManager>(protein, weighted_bins);

        case settings::hist::HistogramManagerChoice::HistogramManagerMT:
            return create_manager<HistogramManagerMT>(protein, weighted_bins);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
            return create_manager<HistogramManagerMTFFAvg>(protein, weighted_bins);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
            return create_manager<HistogramManagerMTFFExplicit>(protein, weighted_bins);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: 
            return std::make_unique<HistogramManagerMTFFGrid>(protein);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface: 
            return std::make_unique<HistogramManagerMTFFGridSurface>(protein);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv: 
            return std::make_unique<HistogramManagerMTFFGridScalableExv>(protein);

        case settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT:
            return create_manager<SymmetryManagerMT>(protein, weighted_bins);

        case settings::hist::HistogramManagerChoice::PartialHistogramManager:
            return create_manager<PartialHistogramManager>(protein, weighted_bins);

        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
            return create_manager<PartialHistogramManagerMT>(protein, weighted_bins);

        case settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT:
            return create_manager<PartialSymmetryManagerMT>(protein, weighted_bins);

        // case settings::hist::HistogramManagerChoice::DebugManager:
        //     return std::make_unique<DebugManager<true>>(protein);

        case settings::hist::HistogramManagerChoice::FoXSManager:
        case settings::hist::HistogramManagerChoice::PepsiManager:
        case settings::hist::HistogramManagerChoice::CrysolManager:
            // FoXSManager, PepsiManager, and CrysolManager are all extensions of the HistogramManagerMTFFExplicit method
            return create_manager<HistogramManagerMTFFExplicit>(protein, weighted_bins);

        default:
            throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
    }
}