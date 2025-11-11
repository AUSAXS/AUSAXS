// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <data/Body.h>
#include <constants/Constants.h>
#include <utility/Concepts.h>

#include <vector>

namespace ausaxs::hist::detail {
    /**
     * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
     *        The idea is that by only extracting the absolute necessities for the distance calculation such as the coordinates and weight,
     *        more values can be stored in the cache at any given time. This is further improved by storing the coordinates as floats instead of doubles.
     *        This is meant as a helper class to DistanceCalculator.
     */
    template<bool variable_bin_width>
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

            std::vector<CompactCoordinatesXYZW>& get_data();
            const std::vector<CompactCoordinatesXYZW>& get_data() const;

            CompactCoordinatesXYZW& operator[](unsigned int i);
            const CompactCoordinatesXYZW& operator[](unsigned int i) const;

        protected: 
            std::vector<CompactCoordinatesXYZW> data;
    };
    static_assert(supports_nothrow_move_v<CompactCoordinates<true>>, "CompactCoordinates should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinates<false>>,"CompactCoordinates should support nothrow move semantics.");
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining

template<bool vbw>
inline ausaxs::hist::detail::CompactCoordinates<vbw>::CompactCoordinates(unsigned int size) : data(size) {}

inline ausaxs::hist::detail::CompactCoordinates::CompactCoordinates(std::vector<Vector3<double>>&& coordinates, double weight) : data(coordinates.size()) {
    std::transform(coordinates.begin(), coordinates.end(), data.begin(), [weight] (const Vector3<double>& v) {return CompactCoordinatesXYZW(v, weight);});
}

template<bool vbw>
inline ausaxs::hist::detail::CompactCoordinates<vbw>::CompactCoordinates(const std::vector<data::AtomFF>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = CompactCoordinatesXYZW(a.coordinates(), a.weight());
    }
}

template<bool vbw>
inline ausaxs::hist::detail::CompactCoordinates<vbw>::CompactCoordinates(const std::vector<data::Body>& bodies) : data(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.size_atom();})) {
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i++] = CompactCoordinatesXYZW(a.coordinates(), a.weight());
        }
    }
}

template<bool vbw>
inline ausaxs::hist::detail::CompactCoordinates<vbw>::CompactCoordinates(const std::vector<data::Water>& atoms) : data(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = CompactCoordinatesXYZW(a.coordinates(), a.weight());
    }
}

template<bool vbw>
inline void ausaxs::hist::detail::CompactCoordinates<vbw>::implicit_excluded_volume(double volume_per_atom) {
    double displaced_charge = constants::charge::density::water*volume_per_atom;
    double charge_per_atom = -displaced_charge;
    std::for_each(data.begin(), data.end(), [charge_per_atom] (CompactCoordinatesXYZW& d) {d.value.w += charge_per_atom;});
}

inline std::vector<ausaxs::hist::detail::CompactCoordinatesXYZW>& ausaxs::hist::detail::CompactCoordinates::get_data() {return data;}

inline const std::vector<ausaxs::hist::detail::CompactCoordinatesXYZW>& ausaxs::hist::detail::CompactCoordinates::get_data() const {return data;}

template<bool vbw>
inline std::size_t ausaxs::hist::detail::CompactCoordinates<vbw>::size() const {return data.size();}

inline ausaxs::hist::detail::CompactCoordinatesXYZW& ausaxs::hist::detail::CompactCoordinates::operator[](unsigned int i) {return data[i];}

inline const ausaxs::hist::detail::CompactCoordinatesXYZW& ausaxs::hist::detail::CompactCoordinates::operator[](unsigned int i) const {return data[i];}