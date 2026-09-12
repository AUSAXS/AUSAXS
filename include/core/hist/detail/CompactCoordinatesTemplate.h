// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <hist/detail/data/CompactCoordinatesXYZFF.h>
#include <data/Body.h>
#include <constants/Constants.h>
#include <utility/Concepts.h>

#include <algorithm>
#include <numeric>
#include <random>
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
     *
     *        By extracting only what the distance calculation needs - the coordinates and a
     *        weight or form factor index - and storing them as floats, far more atoms fit in
     *        cache at once. This is meant as a helper class to DistanceCalculator.
     *
     *        The components are stored as separate arrays rather than interleaved [x, y, z, w]
     *        tuples. Interleaved storage forced every kernel to transpose a block before it could
     *        compute a squared distance, at 8-10 shuffle-port operations per block; with separate
     *        arrays there is nothing to transpose. The class exists solely to feed the kernels, so
     *        the layout is chosen for them - element-wise access is available but is a gather and
     *        belongs only on O(N) paths, never in the pair loop.
     */
    template<CompactCoordinatesType CoordType, bool variable_bin_width>
    class CompactCoordinatesTemplate {
        using Mapper = CoordinateTypeMapper<CoordType>;
        using AtomType = typename Mapper::atom_type;
        using BlockType = typename Mapper::block_type;
        using NonCoordinateType = typename Mapper::non_coordinate_type;

        public:
            CompactCoordinatesTemplate() = default;
            CompactCoordinatesTemplate(const std::vector<data::AtomFF>& body);
            CompactCoordinatesTemplate(const std::vector<data::Body>& bodies);
            CompactCoordinatesTemplate(const std::vector<data::Water>& atoms);

            /**
             * @brief Calculate and subtract the average excluded volume charge from each atom to implicitly account for the excluded volume contribution.
             */
            void implicit_excluded_volume(double volume_per_atom);

            /**
             * @brief Randomly permute the atom order.
             *
             * The distance histogram is a sum over all pairs, so this leaves the result unchanged
             * up to floating-point summation order. Its purpose is to break the spatial
             * correlation that linear file order carries: with atoms stored in file order,
             * consecutive inner-loop atoms sit close to the outer-loop atom, their distances land
             * in the same bin, and the accumulation serialises on store-to-load forwarding. See
             * decorrelate_order() in AtomOrdering.h, which decides when this is worth doing.
             *
             * @param seed Fixed by default so results stay reproducible run to run.
             */
            void shuffle_order(unsigned int seed = 0x9e3779b9u);

            /**
             * @brief Apply @a f to every stored position.
             */
            template<typename F>
            void transform_coordinates(F&& f);

            /**
             * @brief Multiply every stored position by @a scale.
             */
            void scale_coordinates(double scale);

            /**
             * @brief Get the non-coordinate (fourth-component) value.
             */
            NonCoordinateType get_non_coordinate_value(unsigned int i) const;
            NonCoordinateType& get_non_coordinate_value(unsigned int i);

            float x(unsigned int i) const {return _x[i];}
            float y(unsigned int i) const {return _y[i];}
            float z(unsigned int i) const {return _z[i];}

            Vector3<float> position(unsigned int i) const {return {_x[i], _y[i], _z[i]};}
            void set_position(unsigned int i, const Vector3<float>& v) {_x[i] = v.x(); _y[i] = v.y(); _z[i] = v.z();}

            /**
             * @brief The @a c'th component (0=x, 1=y, 2=z, 3=weight/ff) of atom @a i, as a double.
             *        For diagnostics only.
             */
            double component(unsigned int i, unsigned int c) const;

            /**
             * @brief The atom at index @a i, as the kernels take it.
             */
            AtomType atom(unsigned int i) const;

            /**
             * @brief A block of atoms starting at index @a i, as the kernels take it.
             *        The caller must guarantee that the block width it then reads is in bounds.
             */
            BlockType block(unsigned int i) const;

            std::size_t size() const {return _x.size();}
            bool empty() const {return _x.empty();}

        protected:
            std::vector<float> _x, _y, _z;
            std::vector<NonCoordinateType> _w;

        private:
            void resize(std::size_t n);
            void assign(std::size_t i, const data::AtomFF& a);
            void assign(std::size_t i, const data::Water& a);
    };
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::resize(std::size_t n) {
    _x.resize(n);
    _y.resize(n);
    _z.resize(n);
    _w.resize(n);
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::assign(std::size_t i, const data::AtomFF& a) {
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
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::assign(std::size_t i, const data::Water& a) {
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
inline double ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::component(unsigned int i, unsigned int c) const {
    switch (c) {
        case 0: return _x[i];
        case 1: return _y[i];
        case 2: return _z[i];
        default: return static_cast<double>(_w[i]);
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::AtomType
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::atom(unsigned int i) const {
    return AtomType{_x[i], _y[i], _z[i], _w[i]};
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::BlockType
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::block(unsigned int i) const {
    return BlockType{_x.data() + i, _y.data() + i, _z.data() + i, _w.data() + i};
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::NonCoordinateType
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::get_non_coordinate_value(unsigned int i) const {
    return _w[i];
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::NonCoordinateType&
ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::get_non_coordinate_value(unsigned int i) {
    return _w[i];
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(const std::vector<data::AtomFF>& atoms) {
    resize(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {assign(i, atoms[i]);}
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(const std::vector<data::Body>& bodies) {
    resize(std::accumulate(bodies.begin(), bodies.end(), std::size_t{0},
        [] (std::size_t sum, const data::Body& body) {return sum + body.size_atom();}));
    std::size_t i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {assign(i++, a);}
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(const std::vector<data::Water>& atoms) {
    resize(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) {assign(i, atoms[i]);}
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::implicit_excluded_volume(double volume_per_atom) {
    static_assert(std::is_same_v<CoordType, CoordinateTypeXYZW>, "implicit_excluded_volume only works with weight-based CompactCoordinates");
    if constexpr (std::is_same_v<CoordType, CoordinateTypeXYZW>) {
        double displaced_charge = constants::charge::density::water*volume_per_atom;
        float charge_per_atom = static_cast<float>(-displaced_charge);
        std::for_each(_w.begin(), _w.end(), [charge_per_atom] (float& w) {w += charge_per_atom;});
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::shuffle_order(unsigned int seed) {
    // one permutation applied to every component, so an atom stays intact
    const std::size_t n = size();
    std::vector<unsigned int> perm(n);
    std::iota(perm.begin(), perm.end(), 0u);
    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);

    auto permute = [&perm, n] (auto& v) {
        std::decay_t<decltype(v)> out(n);
        for (std::size_t i = 0; i < n; ++i) {out[i] = v[perm[i]];}
        v = std::move(out);
    };
    permute(_x);
    permute(_y);
    permute(_z);
    permute(_w);
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
template<typename F>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::transform_coordinates(F&& f) {
    for (std::size_t i = 0; i < size(); ++i) {
        Vector3<float> v = f(Vector3<float>{_x[i], _y[i], _z[i]});
        _x[i] = v.x();
        _y[i] = v.y();
        _z[i] = v.z();
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::scale_coordinates(double scale) {
    const float f = static_cast<float>(scale);
    for (std::size_t i = 0; i < size(); ++i) {_x[i] *= f; _y[i] *= f; _z[i] *= f;}
}
