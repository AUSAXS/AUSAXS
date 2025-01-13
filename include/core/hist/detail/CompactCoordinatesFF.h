#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesData.h>
#include <form_factor/FormFactorType.h>
#include <utility/Concepts.h>
#include <data/Body.h>

#include <vector>

namespace ausaxs::hist::detail {
    /**
     * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
     *        The idea is that by only extracting the absolute necessities for the distance calculation such as the coordinates and weight,
     *        more values can be stored in the cache at any given time. This is further improved by storing the coordinates as floats instead of doubles.
     *        This is meant as a helper class to DistanceCalculator.
     */
    class CompactCoordinatesFF : public CompactCoordinates {
        public:
            CompactCoordinatesFF() = default;

            /**
             * @brief Extract the necessary coordinates and weights from a body. 
             */
            CompactCoordinatesFF(const data::Body& body);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of bodies. 
             */
            CompactCoordinatesFF(const std::vector<data::Body>& bodies);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
             */
            CompactCoordinatesFF(const std::vector<data::Water>& atoms);

            /**
             * @brief Get the form factor type of the atom at index i.
             */
            unsigned int get_ff_type(unsigned int i) const;

            std::vector<unsigned int> ff_types;
    };
}

//#########################################//
//############ IMPLEMENTATION #############//
//#########################################//

// implementation defined in header to support efficient inlining

inline ausaxs::hist::detail::CompactCoordinatesFF::CompactCoordinatesFF(const data::Body& body) : CompactCoordinates(body.size_atom()), ff_types(body.size_atom()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = body.get_atom(i); 
        data[i] = hist::detail::CompactCoordinatesData(a.coordinates(), a.weight());
        ff_types[i] = static_cast<int>(a.form_factor_type());
    }
}

inline ausaxs::hist::detail::CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::Body>& bodies) 
    : CompactCoordinates(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.size_atom();})) 
{
    ff_types.resize(size());
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            data[i] = hist::detail::CompactCoordinatesData(a.coordinates(), a.weight());
            ff_types[i++] = static_cast<int>(a.form_factor_type());
        }
    }
}

inline ausaxs::hist::detail::CompactCoordinatesFF::CompactCoordinatesFF(const std::vector<data::Water>& atoms) : CompactCoordinates(atoms.size()), ff_types(atoms.size()) {
    for (unsigned int i = 0; i < size(); ++i) {
        const auto& a = atoms[i]; 
        data[i] = hist::detail::CompactCoordinatesData(a.coordinates(), a.weight());
        ff_types[i] = static_cast<int>(a.form_factor_type());
    }
}

inline unsigned int ausaxs::hist::detail::CompactCoordinatesFF::get_ff_type(unsigned int i) const {return ff_types[i];}