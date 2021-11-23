#include "settings.h"
#include <vector>

namespace setting {
    namespace grid {
        PlacementStrategyChoice psc = AxesStrategy; 
        CullingStrategyChoice csc = CounterStrategy;

        double percent_water = 0.1;

        namespace placement {
            double min_score = 0.1; 
        }
    }

    namespace protein {
        double grid_width = 1; 
    }
}