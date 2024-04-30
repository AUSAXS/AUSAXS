/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hydrate/placement/PlacementStrategy.h>
#include <hydrate/GridMember.h>

using namespace grid;

PlacementStrategy::PlacementStrategy(Grid* grid) {this->grid = grid;}

PlacementStrategy::~PlacementStrategy() = default;
