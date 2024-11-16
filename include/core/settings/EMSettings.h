#pragma once

#include <utility/Limit.h>

namespace ausaxs::settings::em {
    extern unsigned int sample_frequency; // How often a bin is sampled in any direction.
    extern double concentration;          // The concentration in mg/mL used when calculating the absolute intensity scale for simulations.
    extern unsigned int charge_levels;    // The number of partial histograms to utilize.

    extern bool hydrate;                  // Whether to hydrate the protein in the EM algorithm.
    extern bool mass_axis;                // Whether to use a mass axis in place of the threshold axis. 

    extern bool save_pdb;                 // Whether to save the final atomic structure as a PDB file.
    extern Limit alpha_levels;            // The range of alpha-levels to search.

    extern bool fixed_weights;            // Whether to use fixed or dynamic weights for the EM algorithm. Fixed weights means that all atoms will have the same weight of 1.
    extern bool plot_landscapes;          // Whether to plot the evaluated chi2 points. Produces 2 plots; one of the full landscape and another of the area near the minimum. The number of points is roughly determined by setting::em::evals

    namespace simulation {
        extern bool noise; // Whether to generate noise for the simulations.
    }
}