#pragma once

#include <hist/detail/CompactCoordinatesData.h>
#include <utility/Concepts.h>
#include <data/DataFwd.h>

#include <vector>

namespace hist::detail {
    /**
     * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
     *        The idea is that by only extracting the absolute necessities for the distance calculation such as the coordinates and weight,
     *        more values can be stored in the cache at any given time. This is further improved by storing the coordinates as floats instead of doubles.
     *        This is meant as a helper class to DistanceCalculator.
     */
    struct CompactCoordinates {
        CompactCoordinates() = default;

        /**
         * @brief Extract the necessary coordinates and weights from a body. 
         */
        CompactCoordinates(const data::Body& body);

        /**
         * @brief Extract the necessary coordinates and weights from a vector of bodies. 
         */
        CompactCoordinates(const std::vector<data::Body>& bodies);

        /**
         * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
         */
        CompactCoordinates(const std::vector<data::record::Water>& atoms);

        unsigned int get_size() const;

        CompactCoordinatesData& operator[](unsigned int i);
        const CompactCoordinatesData& operator[](unsigned int i) const;

        private: 
            unsigned int size;
            std::vector<CompactCoordinatesData> data;
    };
}