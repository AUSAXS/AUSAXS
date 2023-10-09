#include <hydrate/GridMember.h>
#include <data/Water.h>
#include <data/Atom.h>

using namespace grid;

template<typename T>
GridMember<T>::GridMember(const GridMember<T>& gm) : atom(gm.atom), loc(gm.loc), expanded_volume(gm.expanded_volume) {}

template<typename T>
GridMember<T>::GridMember(const GridMember<T>&& gm) noexcept : atom(std::move(gm.atom)), loc(std::move(gm.loc)), expanded_volume(std::move(gm.expanded_volume)) {}

template<typename T>
GridMember<T>::GridMember(const T& atom, Vector3<int> loc) : atom(atom), loc(std::move(loc)) {}

template<typename T>
GridMember<T>::GridMember(T&& atom, Vector3<int> loc) : atom(std::move(atom)), loc(std::move(loc)) {}

template<typename T>
Vector3<int>& GridMember<T>::get_loc() {return loc;}

template<typename T>
const Vector3<int>& GridMember<T>::get_loc() const {return loc;}

template<typename T>
bool GridMember<T>::is_expanded() const {return expanded_volume;}

template<typename T>
void GridMember<T>::set_expanded(bool b) {expanded_volume = b;}

template<typename T>
T& GridMember<T>::get_atom() {return atom;}

template<typename T>
const T& GridMember<T>::get_atom() const {return atom;}

template<typename T>
constants::atom_t GridMember<T>::get_atom_type() const {return atom.get_element();}

template class grid::GridMember<Atom>;
template class grid::GridMember<Water>;