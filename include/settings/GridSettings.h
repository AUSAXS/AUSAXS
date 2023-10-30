#pragma once

#include <utility/Limit3D.h>

namespace settings {
    namespace grid {
        extern double water_scaling; // The number of generated water molecules as a percent of the number of atoms.
        extern double width;         // Width of each bin of the grid used to represent this protein in Å.
        extern double scaling;       // The percent increase in grid size in all dimensions when the grid size is automatically deduced based on an input vector of atoms.
        extern bool cubic;           // Whether to generate a cubic grid. This is primarily intended for rigid body optimization, to ensure there's enough space for all possible conformations.
        extern double rvol;          // The radius of the excluded volume sphere around each atom.

        namespace detail {
            extern double min_score; // (0.5 + min_score) is the minimum percentage of radial lines which must not intersect anything to place a water molecule.
        }
    };
}

namespace settings::grid {
    enum class PlacementStrategy {
        AxesStrategy, 
        RadialStrategy, 
        JanStrategy
    };
    extern PlacementStrategy placement_strategy;
}

namespace settings::grid {
    enum class CullingStrategy {
        CounterStrategy, 
        OutlierStrategy, 
        RandomStrategy
    };
    extern CullingStrategy culling_strategy;
}