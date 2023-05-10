#include <hydrate/culling/RandomCulling.h>

#include <random>

using namespace grid;

RandomCulling::~RandomCulling() = default;

std::vector<Water> RandomCulling::cull(std::vector<GridMember<Water>>& placed_water) const {
    std::shuffle(placed_water.begin(), placed_water.end(), std::mt19937{std::random_device{}()}); // shuffle the molecules
    return CounterCulling::cull(placed_water);
}
