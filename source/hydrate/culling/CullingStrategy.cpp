#include <hydrate/culling/CullingStrategy.h>

using namespace grid;

CullingStrategy::CullingStrategy(Grid* grid) : grid(grid) {}

CullingStrategy::~CullingStrategy() = default;

void CullingStrategy::set_target_count(unsigned int target_count) {this->target_count = target_count;}
