/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/detail/GridObj.h>
#include <math/Vector3.h>

using namespace grid::detail;

GridObj::GridObj(unsigned int x, unsigned int y, unsigned int z) : container::Container3D<State>(x, y, z, EMPTY) {}

State& GridObj::index(const Vector3<int>& v) {return index(v.x(), v.y(), v.z());}
const State& GridObj::index(const Vector3<int>& v) const {return index(v.x(), v.y(), v.z());}

bool GridObj::is_volume(unsigned int x, unsigned int y, unsigned int z) const {return is_volume(index(x, y, z));}
bool GridObj::is_volume(State s) const {return s & VOLUME;}

bool GridObj::is_empty_or_volume(unsigned int x, unsigned int y, unsigned int z) const {return is_empty_or_volume(index(x, y, z));}
bool GridObj::is_empty_or_volume(State s) const {return s & (EMPTY | VOLUME);}

bool GridObj::is_empty(unsigned int x, unsigned int y, unsigned int z) const {return is_empty(index(x, y, z));}
bool GridObj::is_empty(State s) const {return s & EMPTY;}

bool GridObj::is_atom_area(unsigned int x, unsigned int y, unsigned int z) const {return is_atom_area(index(x, y, z));}
bool GridObj::is_atom_area(State s) const {return s & A_AREA;}

bool GridObj::is_atom_area_or_volume(unsigned int x, unsigned int y, unsigned int z) const {return is_atom_area_or_volume(index(x, y, z));}
bool GridObj::is_atom_area_or_volume(State s) const {return s & (A_AREA | VOLUME);}

bool GridObj::is_water_area(unsigned int x, unsigned int y, unsigned int z) const {return is_water_area(index(x, y, z));}
bool GridObj::is_water_area(State s) const {return s & W_AREA;}

bool GridObj::is_atom_center(unsigned int x, unsigned int y, unsigned int z) const {return is_atom_center(index(x, y, z));}
bool GridObj::is_atom_center(State s) const {return s & A_CENTER;}

bool GridObj::is_water_center(unsigned int x, unsigned int y, unsigned int z) const {return is_water_center(index(x, y, z));}
bool GridObj::is_water_center(State s) const {return s & W_CENTER;}