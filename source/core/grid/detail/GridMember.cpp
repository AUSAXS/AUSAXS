/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/detail/GridMember.h>

using namespace ausaxs;
using namespace ausaxs::grid;

template<grid_member_t T>
GridMember<T>::GridMember() = default;

template<grid_member_t T>
GridMember<T>::~GridMember() = default;

template<grid_member_t T>
GridMember<T>::GridMember(const GridMember<T>& gm) : atom(gm.atom), loc(gm.loc), expanded_volume(gm.expanded_volume) {}

template<grid_member_t T>
GridMember<T>::GridMember(const GridMember<T>&& gm) noexcept : atom(std::move(gm.atom)), loc(std::move(gm.loc)), expanded_volume(std::move(gm.expanded_volume)) {}

template<grid_member_t T>
GridMember<T>::GridMember(const T& atom, Vector3<int> loc) : atom(atom), loc(std::move(loc)) {}

template<grid_member_t T>
GridMember<T>::GridMember(T&& atom, Vector3<int> loc) : atom(std::move(atom)), loc(std::move(loc)) {}

template<grid_member_t T>
Vector3<int>& GridMember<T>::get_bin_loc() {return loc;}

template<grid_member_t T>
const Vector3<int>& GridMember<T>::get_bin_loc() const {return loc;}

template<grid_member_t T>
Vector3<double>& GridMember<T>::get_absolute_loc() {return atom.coordinates();}

template<grid_member_t T>
const Vector3<double>& GridMember<T>::get_absolute_loc() const {return atom.coordinates();}

template<grid_member_t T>
bool GridMember<T>::is_expanded() const {return expanded_volume;}

template<grid_member_t T>
void GridMember<T>::set_expanded(bool b) {expanded_volume = b;}

template<grid_member_t T>
T& GridMember<T>::get_atom() {return atom;}

template<grid_member_t T>
const T& GridMember<T>::get_atom() const {return atom;}

template<grid_member_t T>
form_factor::form_factor_t GridMember<T>::get_atom_type() const {return atom.form_factor_type();}

template class grid::GridMember<data::AtomFF>;
template class grid::GridMember<data::Water>;

static_assert(std::is_trivial_v<data::AtomFF>,          "GridMember is not trivial");
static_assert(std::is_standard_layout_v<data::AtomFF>,  "GridMember is not standard layout");
static_assert(supports_nothrow_move_v<data::AtomFF>,    "GridMember should support nothrow move semantics.");