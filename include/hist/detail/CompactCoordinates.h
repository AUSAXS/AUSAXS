#pragma once

#include <utility/Concepts.h>

#include <vector>

template<numeric T> class Vector3;
class Body;
class Water;
namespace hist::detail {
    /**
     * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
     *        The idea is that by only extracting the absolute necessities for the distance calculation such as the coordinates and weight,
     *        more values can be stored in the cache at any given time. This is further improved by storing the coordinates as floats instead of doubles.
     *        This is meant as a helper class to DistanceCalculator.
     */
    struct CompactCoordinates {
        struct Data {
            Data();
            Data(const Vector3<double>& v, float w);
            float x, y, z, w;
        };
        static_assert(sizeof(Data) == 16, "hist::detail::CompactCoordinates::Data is not 16 bytes");

        CompactCoordinates() = default;

        /**
         * @brief Extract the necessary coordinates and weights from a body. 
         */
        CompactCoordinates(const Body& body);

        /**
         * @brief Extract the necessary coordinates and weights from a vector of bodies. 
         */
        CompactCoordinates(const std::vector<Body>& bodies);

        /**
         * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
         */
        CompactCoordinates(const std::vector<Water>& atoms);

        unsigned int get_size() const;

        Data& operator[](unsigned int i);
        const Data& operator[](unsigned int i) const;

        private: 
            unsigned int size;
            std::vector<Data> data;
    };
}