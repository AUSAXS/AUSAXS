// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <math/Transform.h>
#include <utility/Concepts.h>
#include <utility/Random.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <vector>

namespace ausaxs::hist::detail {
    /**
     * @brief A compact representation of the coordinates and weight of all atoms in a body.
     *        This is only designed as a temporary representation for the duration of the histogram calculation.
     */
    class CompactCoordinates {
        public:
            CompactCoordinates() = default;

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
             * @brief Append the contents of @a other.
             */
            void append(const CompactCoordinates& other);

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

            float get_weight(int i) const {return _w[i];}
            float& get_weight(int i) {return _w[i];}

            float x(int i) const {return _x[i];}
            float y(int i) const {return _y[i];}
            float z(int i) const {return _z[i];}

            Vector3<float> position(int i) const {return {_x[i], _y[i], _z[i]};}
            void set_position(int i, const Vector3<float>& v) {_x[i] = v.x(); _y[i] = v.y(); _z[i] = v.z();}

            /**
             * @brief The atom at index @a i, as the kernels take it.
             */
            xyzw::Atom atom(int i) const {return xyzw::Atom{.x=_x[i], .y=_y[i], .z=_z[i], .w=_w[i]};}

            /**
             * @brief A block of atoms starting at index @a i, as the kernels take it.
             *        The caller must guarantee that the block width it then reads is in bounds.
             */
            xyzw::Block block(int i) const {return xyzw::Block{.x=_x.data()+i, .y=_y.data()+i, .z=_z.data()+i, .w=_w.data()+i};}

            int size() const {return static_cast<int>(_x.size());}
            bool empty() const {return _x.empty();}

            /**
             * @brief Resize to @a n entries.
             */
            void resize(int n);

        private:
            std::vector<float> _x, _y, _z, _w;

            template<typename T>
            void assign(int i, const T& a);
    };

    static_assert(supports_nothrow_move_v<CompactCoordinates>, "CompactCoordinates should support nothrow move semantics.");
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining

inline void ausaxs::hist::detail::CompactCoordinates::resize(int n) {
    _x.resize(n);
    _y.resize(n);
    _z.resize(n);
    _w.resize(n);
}

template<typename T>
inline void ausaxs::hist::detail::CompactCoordinates::assign(int i, const T& a) {
    const auto& p = a.coordinates();
    _x[i] = static_cast<float>(p.x());
    _y[i] = static_cast<float>(p.y());
    _z[i] = static_cast<float>(p.z());
    _w[i] = static_cast<float>(a.weight());
}

inline void ausaxs::hist::detail::CompactCoordinates::fill(const std::vector<data::AtomFF>& atoms) {
    resize(static_cast<int>(atoms.size()));
    int i = 0;
    for (const auto& a : atoms) {assign(i++, a);}
}

inline void ausaxs::hist::detail::CompactCoordinates::fill_from_atoms(observer_ptr<const data::Molecule> molecule) {
    resize(molecule->size_atom());
    int i = 0;
    for (const auto& a : molecule->iterate_atoms()) {assign(i++, a);}
}

inline void ausaxs::hist::detail::CompactCoordinates::fill_from_waters(observer_ptr<const data::Molecule> molecule) {
    resize(molecule->size_water());
    int i = 0;
    for (const auto& w : molecule->iterate_waters()) {assign(i++, w);}
}

inline void ausaxs::hist::detail::CompactCoordinates::append(const CompactCoordinates& other) {
    _x.insert(_x.end(), other._x.begin(), other._x.end());
    _y.insert(_y.end(), other._y.begin(), other._y.end());
    _z.insert(_z.end(), other._z.begin(), other._z.end());
    _w.insert(_w.end(), other._w.begin(), other._w.end());
}

inline void ausaxs::hist::detail::CompactCoordinates::implicit_excluded_volume(double volume_per_atom) {
    double displaced_charge = constants::charge::density::water*volume_per_atom;
    auto charge_per_atom = static_cast<float>(-displaced_charge);
    std::ranges::for_each(_w, [charge_per_atom] (float& w) {w += charge_per_atom;});
}

inline void ausaxs::hist::detail::CompactCoordinates::shuffle_order() {
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

inline void ausaxs::hist::detail::CompactCoordinates::transform_coordinates(const ausaxs::transform::Affine& t) {
    for (int i = 0; i < size(); ++i) {
        Vector3<float> v = t({static_cast<double>(_x[i]), static_cast<double>(_y[i]), static_cast<double>(_z[i])});
        _x[i] = v.x();
        _y[i] = v.y();
        _z[i] = v.z();
    }
}

inline void ausaxs::hist::detail::CompactCoordinates::scale_coordinates(double scale) {
    auto f = static_cast<float>(scale);
    for (int i = 0; i < size(); ++i) {_x[i] *= f; _y[i] *= f; _z[i] *= f;}
}
