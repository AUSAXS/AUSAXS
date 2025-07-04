// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <grid/detail/GridObj.h>
#include <math/Vector3.h>

using namespace ausaxs;
using namespace ausaxs::grid::detail;

GridObj::GridObj(int x, int y, int z) : container::Container3D<State>(x, y, z, EMPTY) {}

State& GridObj::index(const Vector3<int>& v) {return index(v.x(), v.y(), v.z());}
const State& GridObj::index(const Vector3<int>& v) const {return index(v.x(), v.y(), v.z());}

bool GridObj::is_empty_or_volume(int x, int y, int z) const {return is_empty_or_volume(index(x, y, z));}
bool GridObj::is_empty_or_volume(State s) const {return s == EMPTY || (s & VOLUME);}

bool GridObj::is_empty_or_volume_or_water(int x, int y, int z) const {return is_empty_or_volume_or_water(index(x, y, z));}
bool GridObj::is_empty_or_volume_or_water(State s) const {return s == EMPTY || (s & (VOLUME | W_AREA | W_CENTER));}

bool GridObj::is_only_empty_or_volume(State s) const {return s == EMPTY || ((s & VOLUME) && !(s & ~VOLUME));}
bool GridObj::is_only_empty_or_volume(int x, int y, int z) const {return is_only_empty_or_volume(index(x, y, z));}

bool GridObj::is_empty_or_water(int x, int y, int z) const {return is_empty_or_water(index(x, y, z));}
bool GridObj::is_empty_or_water(State s) const {return s == EMPTY || (s & (W_AREA | W_CENTER));}

bool GridObj::is_atom_area_or_volume(int x, int y, int z) const {return is_atom_area_or_volume(index(x, y, z));}
bool GridObj::is_atom_area_or_volume(State s) const {return s & (A_AREA | VOLUME);}

bool GridObj::is_only_atom_area_or_volume(int x, int y, int z) const {return is_only_atom_area_or_volume(index(x, y, z));}
bool GridObj::is_only_atom_area_or_volume(State s) const {return (s & (A_AREA | VOLUME)) && !(s & ~(A_AREA | VOLUME));}

bool GridObj::is_volume(int x, int y, int z) const {return is_volume(index(x, y, z));}
bool GridObj::is_volume(State s) const {return s & VOLUME;}

bool GridObj::is_only_volume(int x, int y, int z) const {return is_only_volume(index(x, y, z));}
bool GridObj::is_only_volume(State s) const {return (s & VOLUME) && !(s & ~VOLUME);}

bool GridObj::is_empty(int x, int y, int z) const {return is_empty(index(x, y, z));}
bool GridObj::is_empty(State s) const {return s == EMPTY;}

bool GridObj::is_atom_area(int x, int y, int z) const {return is_atom_area(index(x, y, z));}
bool GridObj::is_atom_area(State s) const {return s & A_AREA;}

bool GridObj::is_water_area(int x, int y, int z) const {return is_water_area(index(x, y, z));}
bool GridObj::is_water_area(State s) const {return s & W_AREA;}

bool GridObj::is_atom_center(int x, int y, int z) const {return is_atom_center(index(x, y, z));}
bool GridObj::is_atom_center(State s) const {return s & A_CENTER;}
bool GridObj::is_only_atom_center(int x, int y, int z) const {return is_only_atom_center(index(x, y, z));}
bool GridObj::is_only_atom_center(State s) const {return (s & A_CENTER) && !(s & ~A_CENTER);}

bool GridObj::is_water_center(int x, int y, int z) const {return is_water_center(index(x, y, z));}
bool GridObj::is_water_center(State s) const {return s & W_CENTER;}