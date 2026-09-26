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
#include <settings/FitSettings.h>
#include <settings/HistogramSettings.h>
#include <settings/InternalState.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::hist::factory;

namespace {
    using Choice = settings::hist::HistogramManagerChoice;
    using ExvMethod = settings::exv::ExvMethod;

    template<template<bool> class MANAGER, typename... Args>
    std::unique_ptr<hist::IHistogramManager> create_manager(bool weighted_bins, observer_ptr<const data::Molecule> protein, Args... args) {
        if (weighted_bins) {
            return std::make_unique<MANAGER<true>>(protein, args...);
        }
        return std::make_unique<MANAGER<false>>(protein, args...);
    }

    /**
     * @brief Whether @a exv_method resolves the atoms by form factor type, rather than weighting each of them.
     */
    bool uses_form_factors(ExvMethod exv_method) {
        switch (exv_method) {
            case ExvMethod::Simple:
            case ExvMethod::None:
                return false;

            // we explicitly write each case to ensure we will get a compiler warning for new models in the future
            case ExvMethod::Average:
            case ExvMethod::Fraser:
            case ExvMethod::Grid:
            case ExvMethod::GridSurface:
            case ExvMethod::GridScalable:
            case ExvMethod::CRYSOL:
            case ExvMethod::FoXS:
            case ExvMethod::Pepsi:
            case ExvMethod::WAXSiS:
                return true;
        }
        throw except::unexpected("hist::factory::uses_form_factors: Unknown ExvMethod. Did you forget to add it to the switch statement?");
    }

    /**
     * @brief The manager of the grid model @a exv_method. These are only implemented for recalculating everything, so the kind
     *        @a choice only decides what to warn about.
     */
    std::unique_ptr<hist::IHistogramManager> create_grid_manager(observer_ptr<const data::Molecule> protein, Choice choice, ExvMethod exv_method) {
        bool partial = choice == Choice::PartialHistogramManager || choice == Choice::PartialHistogramManagerMT || choice == Choice::PartialHistogramSymmetryManagerMT;
        if (partial) {
            console::print_warning(
                "construct_histogram_manager: A partial histogram manager was requested, but the grid excluded volume models have no partial implementation. "
                "Every update will recalculate the full histogram."
            );
        }
        if (protein->symmetry().has_symmetries()) {
            console::print_warning(
                "construct_histogram_manager: Molecule contains symmetries, but the grid excluded volume models do not support them. "
                "Symmetries will be ignored. "
            );
        }

        // without exv fitting, the plain grid manager gives the same result faster
        switch (exv_method) {
            case ExvMethod::GridScalable:
                if (settings::fit::fit_excluded_volume) {return std::make_unique<hist::HistogramManagerMTFFGridScalableExv>(protein);}
                break;
            case ExvMethod::GridSurface:
                if (settings::fit::fit_excluded_volume) {return std::make_unique<hist::HistogramManagerMTFFGridSurface>(protein);}
                break;
            default:
                break;
        }
        return std::make_unique<hist::HistogramManagerMTFFGrid>(protein);
    }
}

bool hist::factory::uses_grid_exv(settings::exv::ExvMethod exv_method) {
    switch (exv_method) {
        case ExvMethod::Grid:
        case ExvMethod::GridSurface:
        case ExvMethod::GridScalable:
        case ExvMethod::WAXSiS:
            return true;
        default:
            return false;
    }
}

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, bool weighted_bins, settings::exv::ExvMethod exv_method
) {
    // symmetry-awareness is derived from the molecule, never from a setting
    auto choice = settings::hist::get_histogram_manager();
    if (protein->symmetry().has_symmetries()) {
        switch (choice) {
            case Choice::HistogramManager:
            case Choice::HistogramManagerMT:
                choice = Choice::HistogramSymmetryManagerMT;
                break;
            case Choice::PartialHistogramManager:
            case Choice::PartialHistogramManagerMT:
                choice = Choice::PartialHistogramSymmetryManagerMT;
                break;
            default:
                break;
        }
    }
    return construct_histogram_manager(protein, choice, weighted_bins, exv_method);
}

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, settings::hist::HistogramManagerChoice choice, bool weighted_bins, settings::exv::ExvMethod exv_method
) {
    if (uses_grid_exv(exv_method)) {return create_grid_manager(protein, choice, exv_method);}

    bool ff = uses_form_factors(exv_method);
    switch (choice) {
        case Choice::HistogramManager:
            if (!ff) {return create_manager<HistogramManager>(weighted_bins, protein);}
            [[fallthrough]]; // the single-threaded reference implementations are weighted only

        case Choice::HistogramManagerMT:
            if (!ff) {return create_manager<HistogramManagerMT>(weighted_bins, protein);}
            if (exv_method == ExvMethod::Average) {return create_manager<HistogramManagerMTFFAvg>(weighted_bins, protein);}
            return create_manager<HistogramManagerMTFFExplicit>(weighted_bins, protein, exv_method);

        case Choice::HistogramSymmetryManagerMT:
            if (!ff) {return create_manager<SymmetryManagerMT>(weighted_bins, protein);}
            return create_manager<SymmetryManagerMTFF>(weighted_bins, protein, exv_method);

        case Choice::PartialHistogramManager:
            if (!ff) {return create_manager<PartialHistogramManager>(weighted_bins, protein);}
            [[fallthrough]]; // the single-threaded reference implementations are weighted only

        case Choice::PartialHistogramManagerMT:
            if (!ff) {return create_manager<PartialHistogramManagerMT>(weighted_bins, protein);}
            return create_manager<PartialHistogramManagerMTFF>(weighted_bins, protein, exv_method);

        case Choice::PartialHistogramSymmetryManagerMT:
            if (!ff) {return create_manager<PartialSymmetryManagerMT>(weighted_bins, protein);}
            return create_manager<PartialSymmetryManagerMTFF>(weighted_bins, protein, exv_method);

        default:
            throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
    }
}