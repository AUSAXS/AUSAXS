#pragma once

namespace settings {
    namespace molecule {
        extern bool center;                 // Decides if the structure will be centered at origo.
        extern bool use_effective_charge;   // Decides whether the charge of the displaced water will be included.
        extern bool throw_on_unknown_atom;  // Decides whether an exception will be thrown if an unknown atom is encountered.
        extern bool implicit_hydrogens;     // Decides whether implicit hydrogens will be added to the structure.
        extern bool use_occupancy;          // Decides whether the occupancy of the atoms will be ignored.

        enum class DisplacedVolumeSet {
            Traube,                         // Traube 1895 as used by CRYSOL, Pepsi-SAXS & FoXS
            Voronoi_implicit_H,             // Voronoi volume with implicit hydrogens
            Voronoi_explicit_H,             // Voronoi volume with explicit hydrogens
            MinimumFluctutation_implicit_H, // Minimum fluctuation volume with implicit hydrogens
            MinimumFluctutation_explicit_H, // Minimum fluctuation volume with explicit hydrogens
            vdw,                            // Volumes calculated from the van der Waals radii

            //? Custom displaced volume set. Make sure to define it first with form_factor::storage::detail::set_custom_displaced_volume_set
            Custom,

            //! Remember to update constants::displaced_volume::standard if this is changed
            Default = Voronoi_implicit_H    // Default displaced volume set
        };
        extern DisplacedVolumeSet displaced_volume_set;
    }

    namespace hydrate {
        enum class HydrationStrategy {
            AxesStrategy, 
            RadialStrategy, 
            JanStrategy,
            PepsiStrategy,
            NoStrategy
        };
        extern HydrationStrategy hydration_strategy;

        enum class CullingStrategy {
            CounterStrategy, 
            BodyCounterStrategy,
            OutlierStrategy, 
            RandomCounterStrategy,
            RandomOutlierStrategy,
            RandomBodyCounterStrategy,
            NoStrategy
        };
        extern CullingStrategy culling_strategy;
    }
}