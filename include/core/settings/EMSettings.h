// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <settings/SettingsHelper.h>
#include <utility/Limit.h>

namespace ausaxs::settings {
    struct EXPORT em {
        static Setting<int> sample_frequency; // How often a bin is sampled in any direction.
        static Setting<double> concentration; // The concentration in mg/mL used when calculating the absolute intensity scale for simulations.
        static Setting<int> charge_levels;    // The number of partial histograms to utilize.

        static Setting<bool> hydrate;         // Whether to hydrate the protein in the EM algorithm.
        static Setting<bool> mass_axis;       // Whether to use a mass axis in place of the threshold axis. 

        static Setting<bool> save_pdb;        // Whether to save the final atomic structure as a PDB file.
        static Setting<Limit> alpha_levels;   // The range of alpha-levels to search.

        static Setting<bool> fixed_weights;   // Whether to use fixed or dynamic weights for the EM algorithm. Fixed weights means that all atoms will have the same weight of 1.
        static Setting<bool> plot_landscapes; // Whether to plot the evaluated chi2 points. Produces 2 plots; one of the full landscape and another of the area near the minimum. The number of points is roughly determined by setting::em::evals

        struct simulation {
            static Setting<bool> noise;       // Whether to generate noise for the simulations.
        };
    };
}