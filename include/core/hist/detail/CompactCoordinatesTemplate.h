// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/detail/data/CompactCoordinatesXYZFF.h>
#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <math/Transform.h>
#include <utility/Concepts.h>
#include <utility/Random.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <numeric>
#include <vector>

namespace ausaxs::hist::detail {
    // Type tags for coordinate types
    struct CoordinateTypeXYZW {};
    struct CoordinateTypeXYZFF {};

    template<typename CoordType>
    concept CompactCoordinatesType =
        std::is_same_v<CoordType, CoordinateTypeXYZW> ||
        std::is_same_v<CoordType, CoordinateTypeXYZFF>;

    // Maps a coordinate type tag onto the kernel namespace's Atom and Block types.
    template<typename CoordType>
    struct CoordinateTypeMapper;

    template<>
    struct CoordinateTypeMapper<CoordinateTypeXYZW> {
        using atom_type = xyzw::Atom;
        using block_type = xyzw::Block;
        using non_coordinate_type = float;
    };

    template<>
    struct CoordinateTypeMapper<CoordinateTypeXYZFF> {
        using atom_type = xyzff::Atom;
        using block_type = xyzff::Block;
        using non_coordinate_type = int32_t;
    };

    /**
     * @brief A compact representation of the coordinates and weight of all atoms in a body.
     *        This is only designed as a temporary representation for the duration of the histogram calculation.  
     */
    template<CompactCoordinatesType CoordType, bool variable_bin_width>
    class CompactCoordinatesTemplate {
        using Mapper = CoordinateTypeMapper<CoordType>;
        using AtomType = typename Mapper::atom_type;
        using BlockType = typename Mapper::block_type;
        using NonCoordinateType = typename Mapper::non_coordinate_type;

        public:
            CompactCoordinatesTemplate() = default;

            /**
             * @brief Replace the contents with @a atoms.
             */
            void fill(const std::vector<data::AtomFF>& atoms);

            /**
             * @brief Replace the contents with the atoms (resp. waters) of @a molecule.
             */
            void fill_from_atoms(observer_ptr<const data::Molecule> molecule);
            void fill_from_waters(observer_ptr<const data::Molecule> molecule);

            /**
             * @brief Calculate and subtract the average excluded volume charge from each atom to implicitly account for the excluded volume contribution.
             */
            void implicit_excluded_volume(double volume_per_atom);

            /**
             * @brief Randomly permute the atom order.
             * 
             * This is designed to break the spatial correlation that linear file order typically carries. Modern CPUs use out-of-order 
             * execution, meaning multiple writes to the same bin can be queued out of order. Breaking the spatial correlation means
             * less contention for the same bins and therefore better parallelization. 
             */
            void shuffle_order();

            /**
             * @brief Apply the rigid transform @a t to every stored position.
             */
            void transform_coordinates(const transform::Affine& t);

            /**
             * @brief Multiply every stored position by @a scale.
             */
            void scale_coordinates(double scale);

            /**
             * @brief Get the non-coordinate (fourth-component) value.
             */
            NonCoordinateType get_non_coordinate_value(int i) const;
            NonCoordinateType& get_non_coordinate_value(int i);

            float x(int i) const {return _x[i];}
            float y(int i) const {return _y[i];}
            float z(int i) const {return _z[i];}

            Vector3<float> position(int i) const {return {_x[i], _y[i], _z[i]};}
            void set_position(int i, const Vector3<float>& v) {_x[i] = v.x(); _y[i] = v.y(); _z[i] = v.z();}

            /**
             * @brief The atom at index @a i, as the kernels take it.
             */
            AtomType atom(int i) const;

            /**
             * @brief A block of atoms starting at index @a i, as the kernels take it.
             *        The caller must guarantee that the block width it then reads is in bounds.
             */
            BlockType block(int i) const;

            int size() const {return static_cast<int>(_x.size());}
            bool empty() const {return _x.empty();}

            /**
             * @brief Resize to @a n entries, keeping the first @a n of whatever is already stored.
             *        Any entries this adds are value-initialized, to be filled through set_position()
             *        and get_non_coordinate_value().
             */
            void resize(int n);

        protected:
            std::vector<float> _x, _y, _z;
            std::vector<NonCoordinateType> _w;

        private:
            void assign(int i, const data::AtomFF& a);
            void assign(int i, const data::Water& a);
    };
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::resize(int n) {
    _x.resize(n);
    _y.resize(n);
    _z.resize(n);
    _w.resize(n);
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::assign(int i, const data::AtomFF& a) {
    const auto& p = a.coordinates();
    _x[i] = static_cast<float>(p.x());
    _y[i] = static_cast<float>(p.y());
    _z[i] = static_cast<float>(p.z());
    if constexpr (std::is_same_v<CoordType, CoordinateTypeXYZW>) {
        _w[i] = static_cast<float>(a.weight());
    } else {
        _w[i] = static_cast<int32_t>(a.form_factor_type());
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::assign(int i, const data::Water& a) {
    const auto& p = a.coordinates();
    _x[i] = static_cast<float>(p.x());
    _y[i] = static_cast<float>(p.y());
    _z[i] = static_cast<float>(p.z());
    if constexpr (std::is_same_v<CoordType, CoordinateTypeXYZW>) {
        _w[i] = static_cast<float>(a.weight());
    } else {
        _w[i] = static_cast<int32_t>(ausaxs::data::Water::form_factor_type());
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::AtomType
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::atom(int i) const {
    return AtomType{_x[i], _y[i], _z[i], _w[i]};
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::BlockType
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::block(int i) const {
    return BlockType{_x.data()+i, _y.data()+i, _z.data()+i, _w.data()+i};
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::NonCoordinateType
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::get_non_coordinate_value(int i) const {
    return _w[i];
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::NonCoordinateType&
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::get_non_coordinate_value(int i) {
    return _w[i];
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::fill(const std::vector<data::AtomFF>& atoms) {
    resize(static_cast<int>(atoms.size()));
    int i = 0;
    for (const auto& a : atoms) {assign(i++, a);}
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::fill_from_atoms(observer_ptr<const data::Molecule> molecule) {
    resize(molecule->size_atom());
    int i = 0;
    for (const auto& a : molecule->iterate_atoms()) {assign(i++, a);}
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::fill_from_waters(observer_ptr<const data::Molecule> molecule) {
    resize(molecule->size_water());
    int i = 0;
    for (const auto& w : molecule->iterate_waters()) {assign(i++, w);}
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::implicit_excluded_volume(double volume_per_atom) {
    static_assert(std::is_same_v<CoordType, CoordinateTypeXYZW>, "implicit_excluded_volume only works with weight-based CompactCoordinates");
    if constexpr (std::is_same_v<CoordType, CoordinateTypeXYZW>) {
        double displaced_charge = constants::charge::density::water*volume_per_atom;
        auto charge_per_atom = static_cast<float>(-displaced_charge);
        std::for_each(_w.begin(), _w.end(), [charge_per_atom] (float& w) {w += charge_per_atom;});
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::shuffle_order() {
    // one permutation applied to every component, so an atom stays intact
    int n = size();
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), random::generator());

    auto permute = [&perm, n] (auto& v) {
        std::decay_t<decltype(v)> out(n);
        for (int i = 0; i < n; ++i) {out[i] = v[perm[i]];}
        v = std::move(out);
    };
    permute(_x);
    permute(_y);
    permute(_z);
    permute(_w);
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::transform_coordinates(const ausaxs::transform::Affine& t) {
    for (int i = 0; i < size(); ++i) {
        Vector3<float> v = t({static_cast<double>(_x[i]), static_cast<double>(_y[i]), static_cast<double>(_z[i])});
        _x[i] = v.x();
        _y[i] = v.y();
        _z[i] = v.z();
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::scale_coordinates(double scale) {
    auto f = static_cast<float>(scale);
    for (int i = 0; i < size(); ++i) {_x[i] *= f; _y[i] *= f; _z[i] *= f;}
}
