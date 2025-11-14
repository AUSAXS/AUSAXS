// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/ExvSettings.h>
#include <settings/FitSettings.h>
#include <settings/Flags.h>

using namespace ausaxs;
using namespace ausaxs::hist::factory;

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, bool weighted_bins, bool variable_bin_width
) {
    auto choice = settings::hist::get_histogram_manager();
    bool has_syms = protein->symmetry().has_symmetries();
    if (has_syms) {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager:
            case settings::hist::HistogramManagerChoice::HistogramManagerMT:
                choice = settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT;
                break;
            default: 
                console::print_warning(
                    "construct_histogram_manager: Molecule contains symmetries, but the chosen excluded volume method does not support them. "
                    "Symmetries will be ignored. "
                );
                break;
        }
    }
    return construct_histogram_manager(protein, choice, weighted_bins, variable_bin_width);
}

template<template<bool, bool> class MANAGER>
std::unique_ptr<hist::IHistogramManager> create_manager(observer_ptr<const data::Molecule> protein, bool weighted_bins, bool variable_bin_width) {
    if (weighted_bins) {
        if (variable_bin_width) {
            return std::make_unique<MANAGER<true, true>>(protein);
        } else {
            return std::make_unique<MANAGER<true, false>>(protein);
        }
    } else {
        if (variable_bin_width) {
            return std::make_unique<MANAGER<false, true>>(protein);
        } else {
            return std::make_unique<MANAGER<false, false>>(protein);
        }
    }
}

template<template<bool> class MANAGER>
std::unique_ptr<hist::IHistogramManager> create_manager(observer_ptr<const data::Molecule> protein, bool variable_bin_width) {
    if (variable_bin_width) {
        return std::make_unique<MANAGER<true>>(protein);
    } else {
        return std::make_unique<MANAGER<false>>(protein);
    }
}

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, settings::hist::HistogramManagerChoice choice, bool weighted_bins, bool variable_bin_width
) {
    switch (choice) {
        case settings::hist::HistogramManagerChoice::HistogramManager:
            return create_manager<HistogramManager>(protein, weighted_bins, variable_bin_width);

        case settings::hist::HistogramManagerChoice::HistogramManagerMT:
            return create_manager<HistogramManagerMT>(protein, weighted_bins, variable_bin_width);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
            return create_manager<HistogramManagerMTFFAvg>(protein, weighted_bins, variable_bin_width);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
            return create_manager<HistogramManagerMTFFExplicit>(protein, weighted_bins, variable_bin_width);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: 
            return create_manager<HistogramManagerMTFFGrid>(protein, variable_bin_width);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface: 
            return create_manager<HistogramManagerMTFFGridSurface>(protein, variable_bin_width);

        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv: 
            return create_manager<HistogramManagerMTFFGridScalableExv>(protein, variable_bin_width);

        case settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT:
            return create_manager<SymmetryManagerMT>(protein, weighted_bins, variable_bin_width);

        case settings::hist::HistogramManagerChoice::PartialHistogramManager:
            return create_manager<PartialHistogramManager>(protein, weighted_bins, variable_bin_width);

        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
            return create_manager<PartialHistogramManagerMT>(protein, weighted_bins, variable_bin_width);

        case settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT:
            return create_manager<PartialSymmetryManagerMT>(protein, weighted_bins, variable_bin_width);

        // case settings::hist::HistogramManagerChoice::DebugManager:
        //     return std::make_unique<DebugManager<true>>(protein);

        case settings::hist::HistogramManagerChoice::FoXSManager:
        case settings::hist::HistogramManagerChoice::PepsiManager:
        case settings::hist::HistogramManagerChoice::CrysolManager:
            // FoXSManager, PepsiManager, and CrysolManager are all extensions of the HistogramManagerMTFFExplicit method
            return create_manager<HistogramManagerMTFFExplicit>(protein, weighted_bins, variable_bin_width);

        default:
            throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
    }
}