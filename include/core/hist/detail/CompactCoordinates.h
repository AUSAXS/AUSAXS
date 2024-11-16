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
            CompactCoordinates(const data::Body& body);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of bodies. 
             */
            CompactCoordinates(const std::vector<data::Body>& bodies);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
             */
            CompactCoordinates(const std::vector<data::record::Water>& atoms);

            std::size_t size() const;

            std::vector<CompactCoordinatesData>& get_data();
            const std::vector<CompactCoordinatesData>& get_data() const;

            CompactCoordinatesData& operator[](unsigned int i);
            const CompactCoordinatesData& operator[](unsigned int i) const;

        protected: 
            std::vector<CompactCoordinatesData> data;
    };
}