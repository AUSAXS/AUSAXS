#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <form_factor/FormFactorType.h>
#include <utility/Concepts.h>
#include <data/DataFwd.h>
#include <math/MathFwd.h>

#include <vector>

namespace hist::detail {
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
            CompactCoordinatesFF(const std::vector<data::record::Water>& atoms);

            /**
             * @brief Get the form factor type of the atom at index i.
             */
            unsigned int get_ff_type(unsigned int i) const;

            std::vector<unsigned int> ff_types;
    };
}