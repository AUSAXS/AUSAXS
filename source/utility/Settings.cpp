#include <utility/Settings.h>
#include <math/Vector3.h>

#include <vector>
#include <fstream>

// Default settings
namespace setting {
    namespace figures {
        std::string format = "pdf";
    }

    namespace grid {
        PlacementStrategyChoice psc = PlacementStrategyChoice::RadialStrategy; 
        CullingStrategyChoice csc = CullingStrategyChoice::CounterStrategy;

        double percent_water = 0.1;
        double ra = 2.4;
        double rh = 1.5;
        double ra_effective = 2.4;
        double width = 1; 
        double scaling = 0.25;
        bool cubic = false;
        Limit3D axes(-250, 250, -250, 250, -250, 250);

        namespace placement {
            double min_score = 0.1; 
        }
    }

    namespace protein {
        bool center = true; 
        bool use_effective_charge = true;
    }

    namespace axes {
        double scattering_intensity_plot_binned_width = 0.1;
        Axis scattering_intensity_plot_axis = {1000, 0.001, 1.001};
    }

    namespace fit {
        double q_low = 0;
        double q_high = 1000;
        unsigned int N = 100; 
    }

    namespace rigidbody {
        TransformationStrategyChoice tsc = RigidTransform;
        ParameterGenerationStrategyChoice pgsc = Simple;
        BodySelectStrategyChoice bssc = RandomSelect;
    }

    namespace em {
        CullingStrategyChoice csc = CullingStrategyChoice::CounterStrategy; 
        unsigned int sample_frequency = 1; 
        double concentration = 2;
        unsigned int charge_levels = 100;

        namespace simulation {
            bool noise = true; // Whether to generate noise for the simulations. 
        }
    }
}

void setting::reader::read(const std::string path) {
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Error in settings::reader::read: Could not open file \"" + path + "\"");}

    std::string line; // placeholder for the current line
    while(getline(input, line)) {
        if (line[0] == '#') {continue;} // # signifies comments
        
    }
}