#pragma once
#include <vector>
#include "math/Vector3.h"

// A small container of the various settings. These should be set *before* their respective classes are instantiated. 
namespace setting {
    namespace grid {
        enum PlacementStrategyChoice {AxesStrategy, RadialStrategy};
        enum CullingStrategyChoice {CounterStrategy, OutlierStrategy};

        extern PlacementStrategyChoice psc; // The choice of placement algorithm.
        extern CullingStrategyChoice csc; // The choice of culling algorithm. 
        extern double percent_water; // The number of generated water molecules as a percent of the number of atoms. 

        extern double ra; // Default radius of protein atoms. 
        extern double rh; // Default radius of water molecules.
        extern double width; // Width of each bin of the grid used to represent this protein.
        extern int bins; // Default number of bins
        extern Vector3 base_point; // Default base point

        namespace placement {
            extern double min_score; // (0.5 + min_score) is the minimum percentage of radial lines which must not intersect anything to place a water molecule
        }
    }

    namespace protein {
    }

    namespace axes {
        const double scattering_intensity_plot_binned_width = 1; // The width of each bin for the scattering plots.
        const std::vector<double> scattering_intensity_plot_axes = {1000, 0.001, 1.001}; // Axes used for the Debye scattering intensity plots.
    }
}