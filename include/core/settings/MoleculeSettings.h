#pragma once

namespace settings {
    namespace molecule {
        extern bool center;                 // Decides if the structure will be centered at origo.
        extern bool use_effective_charge;   // Decides whether the charge of the displaced water will be included.
        extern bool throw_on_unknown_atom;  // Decides whether an exception will be thrown if an unknown atom is encountered.
        extern bool implicit_hydrogens;     // Decides whether implicit hydrogens will be added to the structure.
        extern bool use_occupancy;          // Decides whether the occupancy of the atoms will be ignored.
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