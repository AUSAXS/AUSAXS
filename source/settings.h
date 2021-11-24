#pragma once
#include <vector>

// A small container of the various settings. These should be set *before* their respective classes are instantiated. 
namespace setting {
    namespace grid {
        enum PlacementStrategyChoice {AxesStrategy, RadialStrategy};
        enum CullingStrategyChoice {CounterStrategy, OutlierStrategy};

        extern PlacementStrategyChoice psc; // The choice of placement algorithm.
        extern CullingStrategyChoice csc; // The choice of culling algorithm. 
        extern double percent_water; // The number of generated water molecules as a percent of the number of atoms. 

        const double default_ra = 2; // Default radius of protein atoms. 
        const double default_rh = 1.5; // Default radius of water molecules.

        namespace placement {
            extern double min_score; // (0.5 + min_score) is the minimum percentage of radial lines which must not intersect anything to place a water molecule
        }
    }

    namespace protein {
        extern double grid_width; // Width of each bin of the grid used to represent this protein.
        const std::vector<double> scattering_intensity_plot_axes = {1000, 0.001, 1.001}; // Axes used for the Debye scattering intensity plots.
    }
}