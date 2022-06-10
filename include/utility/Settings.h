#pragma once

#include <vector>
#include <string>

#include <math/Vector3.h>
#include <utility/Axis.h>

// A small container of the various settings. These should be set *before* their respective classes are instantiated. 
namespace setting {
    namespace figures {
        extern std::string format; // The output format.
    }

    namespace grid {
        enum class PlacementStrategyChoice {AxesStrategy, RadialStrategy, JanStrategy};
        enum class CullingStrategyChoice {CounterStrategy, OutlierStrategy};

        extern PlacementStrategyChoice psc; // The choice of placement algorithm.
        extern CullingStrategyChoice csc;   // The choice of culling algorithm. 

        extern double percent_water; // The number of generated water molecules as a percent of the number of atoms. 
        extern double ra;            // Radius of protein atoms. 
        extern double rh;            // Radius of water molecules.
        extern double ra_effective;  // Effective radius of protein atoms. This is based on the volume the average atom effectively occupies. 
        extern double width;         // Width of each bin of the grid used to represent this protein.
        extern double scaling;       // The percent increase in grid size in all dimensions when the grid size is automatically deduced based on an input vector of atoms. 

        extern Limit3D axes;         // Default axes for the grid 

        namespace placement {
            extern double min_score; // (0.5 + min_score) is the minimum percentage of radial lines which must not intersect anything to place a water molecule
        }
    }

    namespace protein {
        extern bool center;               // Decides if the structure will be centered at origo. 
        extern bool use_effective_charge; // Decides whether the charge of the displaced water will be included. 
    }

    namespace axes {
        extern double scattering_intensity_plot_binned_width; // The width of each bin for the scattering plots.
        extern Axis scattering_intensity_plot_axis;           // Axes used for the Debye scattering intensity plots.
    }

    namespace fit {
        extern double q_low;   // Lower limit on the used q-values
        extern double q_high;  // Upper limit on the used q-values
        extern unsigned int N; // Number of points sampled when discretizing a model scattering curve
    }

    namespace rigidbody {
        enum TransformationStrategyChoice {RigidTransform};
        enum ParameterGenerationStrategyChoice {Simple};
        enum BodySelectStrategyChoice {RandomSelect, SequentialSelect};

        extern TransformationStrategyChoice tsc;
        extern ParameterGenerationStrategyChoice pgsc;
        extern BodySelectStrategyChoice bssc;
    }

    namespace em {
        enum class CullingStrategyChoice {NoStrategy, CounterStrategy};
        extern CullingStrategyChoice csc; // The choice of culling algorithm. 

        extern unsigned int sample_frequency;     // How often a bin is sampled in any direction. 
        extern double concentration;              // The concentration in mg/mL used when calculating the absolute intensity scale for simulations.

        namespace simulation {
            extern bool noise; // Whether to generate noise for the simulations. 
        }
    }

    // Simple reader for reading settings from a text file
    class reader {
        static void read(const std::string path);
    };
}