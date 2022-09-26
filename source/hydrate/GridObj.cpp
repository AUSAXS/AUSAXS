#include <hydrate/GridObj.h>

GridObj::GridObj(unsigned int x, unsigned int y, unsigned int z) : xdim(x), ydim(y), zdim(z), grid(x, std::vector<std::vector<T>>(y, std::vector<T>(z, EMPTY))) {}

GridObj::T& GridObj::index(unsigned int x, unsigned int y, unsigned int z) {return grid[x][y][z];}
const GridObj::T& GridObj::index(unsigned int x, unsigned int y, unsigned int z) const {return grid[x][y][z];}

GridObj::T& GridObj::index(const Vector3<int>& v) {return index(v.x(), v.y(), v.z());}
const GridObj::T& GridObj::index(const Vector3<int>& v) const {return index(v.x(), v.y(), v.z());}

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