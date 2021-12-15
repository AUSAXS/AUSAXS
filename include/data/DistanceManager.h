#include "data/Body.h"
#include "data/StateManager.h"
#include "ScatteringHistogram.h"

/**
 * @brief A smart distance calculator which efficiently calculates the scattering histogram.
 */
class DistanceCalculator {
public:

    ScatteringHistogram calc_distance() {}

private:
    vector<bool> changed;
    vector<Body>& bodies;
};