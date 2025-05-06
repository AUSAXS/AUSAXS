/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/expansion/GridExpander.h>
#include <grid/Grid.h>

using namespace ausaxs::grid::expander;

GridExpander::GridExpander(observer_ptr<grid::Grid> grid) : grid(grid) {}
GridExpander::~GridExpander() = default;

double GridExpander::to_x(int i) const {return grid->to_x(i);}
double GridExpander::to_y(int j) const {return grid->to_y(j);}
double GridExpander::to_z(int k) const {return grid->to_z(k);}
int& GridExpander::get_volume() const {return grid->volume;}