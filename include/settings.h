#pragma once
#include <vector>
#include <string>
#include "math/Vector3.h"

using std::string, std::vector;

// A small container of the various settings. These should be set *before* their respective classes are instantiated. 
namespace setting {
    namespace figures {
        extern string format; // The output format.
    }

    namespace grid {
        enum PlacementStrategyChoice {AxesStrategy, RadialStrategy};
        enum CullingStrategyChoice {CounterStrategy, OutlierStrategy};

        extern PlacementStrategyChoice psc; // The choice of placement algorithm.
        extern CullingStrategyChoice csc; // The choice of culling algorithm. 
        extern double percent_water; // The number of generated water molecules as a percent of the number of atoms. 

        extern double ra; // Radius of protein atoms. 
        extern double rh; // Radius of water molecules.
        extern double ra_effective; // Effective radius of protein atoms. This is based on the volume the average atom effectively occupies. 
        extern double width; // Width of each bin of the grid used to represent this protein.
        extern int bins; // Default number of bins
        extern Vector3 base_point; // Default base point

        namespace placement {
            extern double min_score; // (0.5 + min_score) is the minimum percentage of radial lines which must not intersect anything to place a water molecule
        }
    }

    namespace protein {
        extern bool center; // Decides if the structure will be centered at origo. 
    }

    namespace axes {
        extern double scattering_intensity_plot_binned_width; // The width of each bin for the scattering plots.
        extern vector<double> scattering_intensity_plot_axes; // Axes used for the Debye scattering intensity plots.
    }

    namespace fit {
        extern double q_low; // Lower limit on the used q-values
        extern double q_high; // Upper limit on the used q-values
    }

    // Simple reader for reading settings from a text file
    class reader {
        static void read(const string path);
    };
}