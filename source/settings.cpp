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
        PlacementStrategyChoice psc = AxesStrategy; 
        CullingStrategyChoice csc = CounterStrategy;

        double percent_water = 0.1;
        double ra = 2.4;
        double rh = 1.5;
        double width = 1; 
        int bins = 501;
        Vector3 base_point = {-250, -250, -250};

        namespace placement {
            double min_score = 0.1; 
        }
    }

    namespace protein {
    }

    namespace axes {
        double scattering_intensity_plot_binned_width = 1;
        vector<double> scattering_intensity_plot_axes = {1000, 0.001, 1.001};
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