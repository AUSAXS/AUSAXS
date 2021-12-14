#include "settings.h"
#include <vector>
#include "math/Vector3.h"

namespace setting {
    namespace grid {
        PlacementStrategyChoice psc = AxesStrategy; 
        CullingStrategyChoice csc = CounterStrategy;

        double percent_water = 0.1;
        double ra = 1.5;
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
        std::vector<double> scattering_intensity_plot_axes = {1000, 0.001, 1.001};
    }
}