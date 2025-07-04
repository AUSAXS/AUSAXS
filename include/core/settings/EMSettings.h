// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <utility/Limit.h>

namespace ausaxs::settings {
    struct EXPORT em {
        static unsigned int sample_frequency; // How often a bin is sampled in any direction.
        static double concentration;          // The concentration in mg/mL used when calculating the absolute intensity scale for simulations.
        static unsigned int charge_levels;    // The number of partial histograms to utilize.

        static bool hydrate;                  // Whether to hydrate the protein in the EM algorithm.
        static bool mass_axis;                // Whether to use a mass axis in place of the threshold axis. 

        static bool save_pdb;                 // Whether to save the final atomic structure as a PDB file.
        static Limit alpha_levels;            // The range of alpha-levels to search.

        static bool fixed_weights;            // Whether to use fixed or dynamic weights for the EM algorithm. Fixed weights means that all atoms will have the same weight of 1.
        static bool plot_landscapes;          // Whether to plot the evaluated chi2 points. Produces 2 plots; one of the full landscape and another of the area near the minimum. The number of points is roughly determined by setting::em::evals

        struct simulation {
            static bool noise; // Whether to generate noise for the simulations.
        };
    };
}