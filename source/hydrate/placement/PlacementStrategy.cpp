#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/GridMember.h>
#include <data/Water.h>

using namespace grid;

PlacementStrategy::PlacementStrategy(Grid* grid) {this->grid = grid;}

PlacementStrategy::~PlacementStrategy() = default;
