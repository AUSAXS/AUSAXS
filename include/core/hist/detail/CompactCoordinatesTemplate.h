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
    template<typename C>
    concept CompactCoordinatesDataType = 
        std::is_same_v<C, xyzw::CompactCoordinatesXYZW> ||
        std::is_same_v<C, xyzff::CompactCoordinatesXYZFF>
    ;

    /**
     * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
     *        The idea is that by only extracting the absolute necessities for the distance calculation such as the coordinates and weight,
     *        more values can be stored in the cache at any given time. This is further improved by storing the coordinates as floats instead of doubles.
     *        This is meant as a helper class to DistanceCalculator.
     */
    template<CompactCoordinatesDataType T>
    class CompactCoordinates {
        public:
            CompactCoordinates() = default;

            /**
             * @brief Create a compact coordinate representation of the given coordinates.
             *        The weights are assumed to be unity.
             */
            CompactCoordinates(std::vector<Vector3<double>>&& coordinates, double weight);

            /**
             * @brief Create a CompactCoordinates object with a given size.
             */
            CompactCoordinates(unsigned int size);

            /**
             * @brief Extract the necessary coordinates and weights from a body. 
             */
            CompactCoordinates(const std::vector<data::Atom>& body);
            CompactCoordinates(const std::vector<data::AtomFF>& body);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of bodies. 
             */
            CompactCoordinates(const std::vector<data::Body>& bodies);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
             */
            CompactCoordinates(const std::vector<data::Water>& atoms);

            /**
             * @brief Calculate and subtract the average excluded volume charge from each atom to implicitly account for the excluded volume contribution.
             */
            void implicit_excluded_volume(double volume_per_atom);

            std::size_t size() const;

            std::vector<T>& get_data();
            const std::vector<T>& get_data() const;

            T& operator[](unsigned int i);
            const T& operator[](unsigned int i) const;

        protected: 
            std::vector<T> data;

        static_assert(supports_nothrow_move_v<CompactCoordinates>, "CompactCoordinates should support nothrow move semantics.");
    };
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining
template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline ausaxs::hist::detail::CompactCoordinates<T>::CompactCoordinates(unsigned int size) : data(size) {}

namespace ausaxs::hist::detail {
    template<ausaxs::hist::detail::CompactCoordinatesDataType T>
    T GenericConstructor(const data::AtomFF& a) {
        if constexpr (std::is_same_v<T, xyzw::CompactCoordinatesXYZW>) {
            return xyzw::CompactCoordinatesXYZW(a.coordinates(), a.weight());
        } else {
            static_assert(std::is_same_v<T, xyzff::CompactCoordinatesXYZFF>, "Type T must be CompactCoordinatesXYZFF");
            return xyzff::CompactCoordinatesXYZFF(a.coordinates(), static_cast<int32_t>(a.form_factor_type()));
        }
    }
}


template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline ausaxs::hist::detail::CompactCoordinates<T>::CompactCoordinates(std::vector<Vector3<double>>&& coordinates, double weight) : data(coordinates.size()) {
    std::transform(coordinates.begin(), coordinates.end(), data.begin(), [weight] (const Vector3<double>& v) {return CompactCoordinatesXYZW(v, weight);});
}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline ausaxs::hist::detail::CompactCoordinates<T>::CompactCoordinates(const std::vector<data::AtomFF>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = GenericConstructor<T>(a);
    }
}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline ausaxs::hist::detail::CompactCoordinates<T>::CompactCoordinates(const std::vector<data::Body>& bodies) : data(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.size_atom();})) {
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i++] = GenericConstructor<T>(a);
        }
    }
}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline ausaxs::hist::detail::CompactCoordinates<T>::CompactCoordinates(const std::vector<data::Water>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = GenericConstructor<T>(a);
    }
}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline void ausaxs::hist::detail::CompactCoordinates<T>::implicit_excluded_volume(double volume_per_atom) {
    double displaced_charge = constants::charge::density::water*volume_per_atom;
    double charge_per_atom = -displaced_charge;
    std::for_each(data.begin(), data.end(), [charge_per_atom] (CompactCoordinatesXYZW& d) {d.value.w += charge_per_atom;});
}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline std::vector<ausaxs::hist::detail::CompactCoordinatesXYZW>& ausaxs::hist::detail::CompactCoordinates<T>::get_data() {return data;}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline const std::vector<ausaxs::hist::detail::CompactCoordinatesXYZW>& ausaxs::hist::detail::CompactCoordinates<T>::get_data() const {return data;}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline std::size_t ausaxs::hist::detail::CompactCoordinates<T>::size() const {return data.size();}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline ausaxs::hist::detail::CompactCoordinatesXYZW& ausaxs::hist::detail::CompactCoordinates<T>::operator[](unsigned int i) {return data[i];}

template<ausaxs::hist::detail::CompactCoordinatesDataType T>
inline const ausaxs::hist::detail::CompactCoordinatesXYZW& ausaxs::hist::detail::CompactCoordinates<T>::operator[](unsigned int i) const {return data[i];}