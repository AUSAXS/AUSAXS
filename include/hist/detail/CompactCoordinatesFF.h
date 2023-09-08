#pragma once

#include <utility/Concepts.h>

#include <vector>
#include <cstdint>

template<numeric T> class Vector3;
class Body;
class Water;
namespace hist {
    namespace detail {
        // The form factor type of an atom. This is intended to be used as an index for best performance.
        enum class form_factor_type : std::uint32_t {
            HYDROGEN,
            CARBON,
            NITROGEN,
            OXYGEN,
            OTHER
        };

        /**
         * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
         *        The idea is that by only extracting the absolute necessities for the distance calculation, more values can be stored
         *        in the cache at any given time. This is meant as a helper class to DistanceCalculator.
         */
        struct CompactCoordinatesFF {
            struct Data {
                Data();
                Data(const Vector3<double>& v, float w);
                float x, y, z, w; 
                form_factor_type type;
            };
            static_assert(sizeof(Data) == 20, "hist::detail::CompactCoordinatesFF::Data is not 17 bytes");

            CompactCoordinatesFF() = default;

            /**
             * @brief Extract the necessary coordinates and weights from a body. 
             */
            CompactCoordinatesFF(const Body& body);

            /**
             * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
             */
            CompactCoordinatesFF(const std::vector<Water>& atoms);

            unsigned int size;
            std::vector<Data> data;
        };
    }
}