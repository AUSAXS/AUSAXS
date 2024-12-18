/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <grid/detail/GridMember.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>

using namespace ausaxs;
using namespace ausaxs::grid;

template<typename T>
GridMember<T>::GridMember(const GridMember<T>& gm) : atom(gm.atom), loc(gm.loc), expanded_volume(gm.expanded_volume) {}

template<typename T>
GridMember<T>::GridMember(const GridMember<T>&& gm) noexcept : atom(std::move(gm.atom)), loc(std::move(gm.loc)), expanded_volume(std::move(gm.expanded_volume)) {}

template<typename T>
GridMember<T>::GridMember(const T& atom, Vector3<int> loc) : atom(atom), loc(std::move(loc)) {}

template<typename T>
GridMember<T>::GridMember(T&& atom, Vector3<int> loc) : atom(std::move(atom)), loc(std::move(loc)) {}

template<typename T>
Vector3<int>& GridMember<T>::get_bin_loc() {return loc;}

template<typename T>
const Vector3<int>& GridMember<T>::get_bin_loc() const {return loc;}

template<typename T>
Vector3<double>& GridMember<T>::get_absolute_loc() {return atom.coordinates();}

template<typename T>
const Vector3<double>& GridMember<T>::get_absolute_loc() const {return atom.coordinates();}

template<typename T>
bool GridMember<T>::is_expanded() const {return expanded_volume;}

template<typename T>
void GridMember<T>::set_expanded(bool b) {expanded_volume = b;}

template<typename T>
T& GridMember<T>::get_atom() {return atom;}

template<typename T>
const T& GridMember<T>::get_atom() const {return atom;}

template<typename T>
form_factor::form_factor_t GridMember<T>::get_atom_type() const {return atom.form_factor_type();}

template class grid::GridMember<data::AtomFF>;
template class grid::GridMember<data::Water>;

static_assert(std::is_trivial_v<data::AtomFF>,          "GridMember is not trivial");
static_assert(std::is_standard_layout_v<data::AtomFF>,  "GridMember is not standard layout");
static_assert(supports_nothrow_move_v<data::AtomFF>,    "GridMember should support nothrow move semantics.");