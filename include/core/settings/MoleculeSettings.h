// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <settings/ExportMacro.h>

namespace ausaxs::settings {
    struct EXPORT molecule {
        static bool center;                 // Decides if the structure will be centered at origo.
        static bool throw_on_unknown_atom;  // Decides whether an exception will be thrown if an unknown atom is encountered.
        static bool allow_unknown_residues; // Decides whether to allow unknown residues (WARNING: implicit hydrogens will be skipped for these residues).
        static bool implicit_hydrogens;     // Decides whether implicit hydrogens will be added to the structure.
        static bool use_occupancy;          // Decides whether the occupancy of the atoms will be ignored.

        enum class ExvSet {
            Traube,                         // Traube 1895 as used by CRYSOL, Pepsi-SAXS & FoXS
            Voronoi_implicit_H,             // Voronoi volume with implicit hydrogens from Schaefer et al.
            Voronoi_explicit_H,             // Voronoi volume with explicit hydrogens from Schaefer et al.
            MinimumFluctutation_implicit_H, // Minimum fluctuation volume with implicit hydrogens from Schaefer et al.
            MinimumFluctutation_explicit_H, // Minimum fluctuation volume with explicit hydrogens from Schaefer et al.
            vdw,                            // Volumes calculated from the van der Waals radii

            //? Custom displaced volume set. Make sure to define it first with form_factor::storage::detail::set_custom_displaced_volume_set
            Custom,

            //! Remember to update constants::displaced_volume::standard if this is changed
            Default = MinimumFluctutation_implicit_H // Default displaced volume set
        };
        static ExvSet exv_set;
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
        static HydrationStrategy hydration_strategy;

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
        static CullingStrategy culling_strategy;

        // Correction added to the sum of van der Waals radii when calculating the distance to the placed water molecules.
        // By default this is tuned to roughly match MD density profiles. 
        static double shell_correction;
    };
}