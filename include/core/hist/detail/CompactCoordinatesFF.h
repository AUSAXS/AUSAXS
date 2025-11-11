// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/data/CompactCoordinatesXYZFF.h>
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
    template<bool variable_bin_width>
    class CompactCoordinatesFF : public CompactCoordinates<variable_bin_width> {
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

template<bool vbw>
void validate_ff_info(const ausaxs::hist::detail::CompactCoordinatesFF<vbw>& data) {
    for (unsigned int i = 0; i < data.size(); ++i) {
        unsigned int ff_type = data.ff_types[i];
        if (ff_type == static_cast<unsigned int>(ausaxs::form_factor::form_factor_t::UNKNOWN)) {
            throw std::runtime_error(
                "CompactCoordinatesFF: Attempted to use an atom with UNKNOWN form factor type.\n"
                "Form factor information is required for the selected excluded volume model."
            );
        }
    }
}

template<bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesFF<vbw>::CompactCoordinatesFF(const data::Body& body) : CompactCoordinates<vbw>(body.size_atom()), ff_types(body.size_atom()) {
    for (unsigned int i = 0; i < this->size(); ++i) {
        const auto& a = body.get_atom(i); 
        this->data[i] = hist::detail::CompactCoordinatesData<vbw>(a.coordinates(), a.weight());
        ff_types[i] = static_cast<int>(a.form_factor_type());
    }
    validate_ff_info(*this);
}

template<bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesFF<vbw>::CompactCoordinatesFF(const std::vector<data::Body>& bodies) 
    : CompactCoordinates<vbw>(std::accumulate(bodies.begin(), bodies.end(), 0, [](unsigned int sum, const data::Body& body) {return sum + body.size_atom();})) 
{
    ff_types.resize(this->size());
    unsigned int i = 0;
    for (const auto& body : bodies) {
        for (const auto& a : body.get_atoms()) {
            this->data[i] = hist::detail::CompactCoordinatesData<vbw>(a.coordinates(), a.weight());
            ff_types[i++] = static_cast<int>(a.form_factor_type());
        }
    }
    validate_ff_info(*this);
}

template<bool vbw>
inline ausaxs::hist::detail::CompactCoordinatesFF<vbw>::CompactCoordinatesFF(const std::vector<data::Water>& atoms) : CompactCoordinates<vbw>(atoms.size()), ff_types(atoms.size()) {
    for (unsigned int i = 0; i < this->size(); ++i) {
        const auto& a = atoms[i]; 
        this->data[i] = hist::detail::CompactCoordinatesData<vbw>(a.coordinates(), a.weight());
        ff_types[i] = static_cast<int>(a.form_factor_type());
    }
    validate_ff_info(*this);
}

template<bool vbw>
inline unsigned int ausaxs::hist::detail::CompactCoordinatesFF<vbw>::get_ff_type(unsigned int i) const {return ff_types[i];}