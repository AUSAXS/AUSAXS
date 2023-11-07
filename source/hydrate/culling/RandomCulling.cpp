#include <hydrate/culling/RandomCulling.h>
#include <hydrate/GridMember.h>
#include <data/record/Water.h>

#include <random>

using namespace grid;
using namespace data::record;

RandomCulling::~RandomCulling() = default;

std::vector<data::record::Water> RandomCulling::cull(std::vector<GridMember<Water>>& placed_water) const {
    auto rng = std::mt19937{std::random_device{}()};
    std::shuffle(placed_water.begin(), placed_water.end(), rng); // shuffle the molecules
    return CounterCulling::cull(placed_water);
}
