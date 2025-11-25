// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <hist/detail/data/CompactCoordinatesXYZFF.h>
#include <data/Body.h>
#include <constants/Constants.h>
#include <utility/Concepts.h>

#include <vector>

namespace ausaxs::hist::detail {
    // Type tags for coordinate types
    struct CoordinateTypeXYZW {};
    struct CoordinateTypeXYZFF {};

    template<typename CoordType>
    concept CompactCoordinatesType = 
        std::is_same_v<CoordType, CoordinateTypeXYZW> ||
        std::is_same_v<CoordType, CoordinateTypeXYZFF>;

    // Helper to get the actual coordinate type
    template<typename CoordType, bool variable_bin_width>
    struct CoordinateTypeMapper;

    template<bool vbw>
    struct CoordinateTypeMapper<CoordinateTypeXYZW, vbw> {
        using type = CompactCoordinatesXYZW<vbw>;
    };

    template<bool vbw>
    struct CoordinateTypeMapper<CoordinateTypeXYZFF, vbw> {
        using type = CompactCoordinatesXYZFF<vbw>;
    };

    /**
     * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
     *        The idea is that by only extracting the absolute necessities for the distance calculation such as the coordinates and weight,
     *        more values can be stored in the cache at any given time. This is further improved by storing the coordinates as floats instead of doubles.
     *        This is meant as a helper class to DistanceCalculator.
     */
    template<CompactCoordinatesType CoordType, bool variable_bin_width>
    class CompactCoordinatesTemplate {
        using DataType = typename CoordinateTypeMapper<CoordType, variable_bin_width>::type;
        using NonCoordinateType = typename std::conditional_t<std::is_same_v<CoordType, CoordinateTypeXYZW>,
            float,
            int32_t
        >;

        public:
            CompactCoordinatesTemplate() = default;

            /**
             * @brief Create a compact coordinate representation of the given coordinates.
             *        The weights are assumed to be unity.
             */
            CompactCoordinatesTemplate(std::vector<Vector3<double>>&& coordinates, double weight);

            /**
             * @brief Create a CompactCoordinates object with a given size.
             */
            CompactCoordinatesTemplate(unsigned int size);

            /**
             * @brief Extract the necessary coordinates and weights from a body. 
             */
            CompactCoordinatesTemplate(const std::vector<data::Atom>& body);
            CompactCoordinatesTemplate(const std::vector<data::AtomFF>& body);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of bodies. 
             */
            CompactCoordinatesTemplate(const std::vector<data::Body>& bodies);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
             */
            CompactCoordinatesTemplate(const std::vector<data::Water>& atoms);

            /**
             * @brief Calculate and subtract the average excluded volume charge from each atom to implicitly account for the excluded volume contribution.
             */
            void implicit_excluded_volume(double volume_per_atom);

            /**
             * @brief Get the non-coordinate (fourth-column) value. 
             */
            NonCoordinateType get_non_coordinate_value(unsigned int i) const;

            std::size_t size() const;

            std::vector<DataType>& get_data();
            const std::vector<DataType>& get_data() const;

            DataType& operator[](unsigned int i);
            const DataType& operator[](unsigned int i) const;

        protected: 
            std::vector<DataType> data;
    };
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining
namespace ausaxs::hist::detail {
    template<CompactCoordinatesType CoordType, bool vbw>
    typename CoordinateTypeMapper<CoordType, vbw>::type GenericConstructor(const data::AtomFF& a) {
        using DataType = typename CoordinateTypeMapper<CoordType, vbw>::type;
        if constexpr (std::is_same_v<CoordType, CoordinateTypeXYZW>) {
            return DataType(a.coordinates(), a.weight());
        } else {
            static_assert(std::is_same_v<CoordType, CoordinateTypeXYZFF>, "Type must be CoordinateTypeXYZFF");
            return DataType(a.coordinates(), static_cast<int32_t>(a.form_factor_type()));
        }
    }

    template<CompactCoordinatesType CoordType, bool vbw>
    typename CoordinateTypeMapper<CoordType, vbw>::type GenericConstructor(const data::Water& a) {
        using DataType = typename CoordinateTypeMapper<CoordType, vbw>::type;
        if constexpr (std::is_same_v<CoordType, CoordinateTypeXYZW>) {
            return DataType(a.coordinates(), a.weight());
        } else {
            static_assert(std::is_same_v<CoordType, CoordinateTypeXYZFF>, "Type must be CoordinateTypeXYZFF");
            return DataType(a.coordinates(), static_cast<int32_t>(a.form_factor_type()));
        }
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(unsigned int size) : data(size) {}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::NonCoordinateType ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::get_non_coordinate_value(unsigned int i) const {
    return data[i].data[3];
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(std::vector<Vector3<double>>&& coordinates, double weight) : data(coordinates.size()) {
    std::transform(coordinates.begin(), coordinates.end(), data.begin(), [weight] (const Vector3<double>& v) {return DataType(v, weight);});
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(const std::vector<data::AtomFF>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = GenericConstructor<CoordType, vbw>(a);
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(const std::vector<data::Body>& bodies) : data(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.size_atom();})) {
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i++] = GenericConstructor<CoordType, vbw>(a);
        }
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::CompactCoordinatesTemplate(const std::vector<data::Water>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = GenericConstructor<CoordType, vbw>(a);
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline void ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::implicit_excluded_volume(double volume_per_atom) {
    static_assert(std::is_same_v<CoordType, CoordinateTypeXYZW>, "implicit_excluded_volume only works with weight-based CompactCoordinates");
    if constexpr (std::is_same_v<CoordType, CoordinateTypeXYZW>) {
        double displaced_charge = constants::charge::density::water*volume_per_atom;
        double charge_per_atom = -displaced_charge;
        std::for_each(data.begin(), data.end(), [charge_per_atom] (DataType& d) {d.value.w += charge_per_atom;});
    }
}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline std::vector<typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::DataType>& ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::get_data() {return data;}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline const std::vector<typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::DataType>& ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::get_data() const {return data;}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline std::size_t ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::size() const {return data.size();}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::DataType& ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::operator[](unsigned int i) {return data[i];}

template<ausaxs::hist::detail::CompactCoordinatesType CoordType, bool vbw>
inline const typename ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::DataType& ausaxs::hist::detail::CompactCoordinatesTemplate<CoordType, vbw>::operator[](unsigned int i) const {return data[i];}