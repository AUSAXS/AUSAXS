// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>
#include <settings/SettingsHelper.h>

namespace ausaxs::settings {
    struct EXPORT molecule {
        static Setting<bool> allow_unknown_atoms;    // Decides whether an exception will be thrown if an unknown atom is encountered.
        static Setting<bool> allow_unknown_residues; // Decides whether to allow unknown residues (WARNING: implicit hydrogens will be skipped for these residues).
        static Setting<bool> implicit_hydrogens;     // Decides whether implicit hydrogens will be added to the structure.
        static Setting<bool> use_occupancy;          // Decides whether the occupancy of the atoms will be ignored.
    };

    struct EXPORT hydrate {
        enum class HydrationStrategy {
            // This strategy attempts to place water molecules along each axis of every atom
            AxesStrategy, 

            // This strategy is a generalization of the AxesStrategy, where water molecules are placed along radial lines from each atom
            RadialStrategy, 

            // This is a simple strategy where the entire grid is scanned for empty positions
            JanStrategy,

            // This is a mimic of the Pepsi hydration method
            PepsiStrategy,

            // No hydration molecules will be placed with this strategy
            NoStrategy
        };
        static Setting<HydrationStrategy> hydration_strategy;

        enum class CullingStrategy {
            // Every N molecule is removed from the list of possible waters to achieve the target number of molecules
            CounterStrategy, 

            // Molecules are removed from the list of possible waters in an attempt to make a more uniform distribution
            OutlierStrategy, 

            // Identical to the CounterStrategy, except the list of possible waters is randomized
            RandomCounterStrategy,

            // Identical to the BodyCounterStrategy, except the list of possible waters is randomized
            RandomOutlierStrategy,

            // No molecules will be removed with this strategy
            NoStrategy
        };
        static Setting<CullingStrategy> culling_strategy;

        // Correction added to the sum of van der Waals radii when calculating the distance to the placed water molecules.
        // By default this is tuned to roughly match MD density profiles. 
        static Setting<double> shell_correction;
    };
}