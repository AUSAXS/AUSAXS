#include "settings.h"
#include "math/Vector3.h"

#include <vector>
#include <fstream>

using std::string, std::vector;

// Default settings
namespace setting {
    namespace figures {
        string format = "pdf";
    }

    namespace grid {
        PlacementStrategyChoice psc = PlacementStrategyChoice::AxesStrategy; 
        CullingStrategyChoice csc = CullingStrategyChoice::CounterStrategy;

        double percent_water = 0.1;
        double ra = 2.4;
        double rh = 1.5;
        double ra_effective = 2.4;
        double width = 1; 
        double scaling = 0.25;
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
        double q_low = 0; // lower limit on the q value
        double q_high = 1000; // upper limit on the q value
    }

    namespace rigidbody {
        TransformationStrategyChoice tsc = RigidTransform;
        ParameterGenerationStrategyChoice pgsc = Simple;
        BodySelectStrategyChoice bssc = RandomSelect;
    }

    namespace em {
        unsigned int max_atoms = 50000;
    }
}

void setting::reader::read(const string path) {
    std::ifstream input(path);
    if (!input.is_open()) {throw std::ios_base::failure("Error in settings::reader::read: Could not open file \"" + path + "\"");}

    string line; // placeholder for the current line
    while(getline(input, line)) {
        if (line[0] == '#') {continue;} // # signifies comments
        
    }
}