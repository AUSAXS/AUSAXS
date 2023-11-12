#include <hydrate/GridObj.h>
#include <math/Vector3.h>

using namespace grid::detail;

GridObj::GridObj(unsigned int x, unsigned int y, unsigned int z) : container::Container3D<State>(x, y, z, EMPTY) {}

// GridObj::State& GridObj::index(unsigned int x, unsigned int y, unsigned int z) {return grid.index(x, y, z);}
// const GridObj::State& GridObj::index(unsigned int x, unsigned int y, unsigned int z) const {return grid.index(x, y, z);}

State& GridObj::index(const Vector3<int>& v) {return index(v.x(), v.y(), v.z());}
const State& GridObj::index(const Vector3<int>& v) const {return index(v.x(), v.y(), v.z());}

bool GridObj::is_volume(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & VOLUME;}
bool GridObj::is_volume(State s) const {return s & VOLUME;}
bool GridObj::is_empty_or_volume(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & (EMPTY | VOLUME);}
bool GridObj::is_empty_or_volume(State s) const {return s & (EMPTY | VOLUME);}
bool GridObj::is_empty(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & EMPTY;}
bool GridObj::is_empty(State s) const {return s & EMPTY;}
bool GridObj::is_atom_area(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & A_AREA;}
bool GridObj::is_atom_area(State s) const {return s & A_AREA;}
bool GridObj::is_atom_area_or_volume(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & (A_AREA | VOLUME);}
bool GridObj::is_atom_area_or_volume(State s) const {return s & (A_AREA | VOLUME);}
bool GridObj::is_water_area(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & W_AREA;}
bool GridObj::is_water_area(State s) const {return s & W_AREA;}
bool GridObj::is_atom_center(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & A_CENTER;}
bool GridObj::is_atom_center(State s) const {return s & A_CENTER;}
bool GridObj::is_water_center(unsigned int x, unsigned int y, unsigned int z) const {return index(x, y, z) & W_CENTER;}
bool GridObj::is_water_center(State s) const {return s & W_CENTER;}

// GridObj::T& GridObj::index(const Vector3<int>& v) {
//     if (v.x() < 0 || v.y() < 0 || v.z() < 0 || v.x() >= xdim || v.y() >= ydim || v.z() >= zdim)
//         throw std::out_of_range("GridObj::index(const Vector3<int>& v) const: Index out of range. (" + std::to_string(v.x()) + ", " + std::to_string(v.y()) + ", " + std::to_string(v.z()) + ")");
//     return index(v.x(), v.y(), v.z());
// }

// const GridObj::T& GridObj::index(const Vector3<int>& v) const {
//     if (v.x() < 0 || v.y() < 0 || v.z() < 0 || v.x() >= xdim || v.y() >= ydim || v.z() >= zdim)
//         throw std::out_of_range("GridObj::index(const Vector3<int>& v) const: Index out of range. (" + std::to_string(v.x()) + ", " + std::to_string(v.y()) + ", " + std::to_string(v.z()) + ")");
//     return index(v.x(), v.y(), v.z());
// }