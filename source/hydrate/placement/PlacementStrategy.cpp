#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/GridMember.h>

using namespace grid;

PlacementStrategy::PlacementStrategy(Grid* grid) {this->grid = grid;}

PlacementStrategy::~PlacementStrategy() = default;
