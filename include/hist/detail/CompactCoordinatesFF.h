#pragma once

#include <hist/detail/FormFactorType.h>
#include <utility/Concepts.h>

#include <vector>

template<numeric T> class Vector3;
class Body;
class Water;
namespace hist::detail {
    /**
     * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
     *        The idea is that by only extracting the absolute necessities for the distance calculation, more values can be stored
     *        in the cache at any given time. This is meant as a helper class to DistanceCalculator.
     */
    struct CompactCoordinatesFF {
        struct Data {
            Data();
            Data(const Vector3<double>& v, float w, hist::detail::ff::form_factor_t ff_type);
            float x, y, z, w; 
            unsigned int ff_type;
        };
        static_assert(sizeof(Data) == 20, "hist::detail::CompactCoordinatesFF::Data is not 20 bytes");

        CompactCoordinatesFF() = default;

        /**
         * @brief Extract the necessary coordinates and weights from a body. 
         */
        CompactCoordinatesFF(const Body& body);

        /**
         * @brief Extract the necessary coordinates and weights from a vector of bodies. 
         */
        CompactCoordinatesFF(const std::vector<Body>& bodies);

        /**
         * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
         */
        CompactCoordinatesFF(const std::vector<Water>& atoms);

        unsigned int size;
        std::vector<Data> data;
    };
}