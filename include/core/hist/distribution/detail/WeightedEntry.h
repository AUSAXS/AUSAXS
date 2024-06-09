#pragma once

#include <constants/Axes.h>
#include <iosfwd>

namespace hist {
    namespace detail {
        /**
         * @brief This struct is a small wrapper around a binned value, automatically keeping tracking the weighted center of the bin. 
         *        This is used to keep the value and weight tracking variables close together in memory, thus improving cache locality for the distance calculations.
         */
        struct WeightedEntry {
            WeightedEntry();
            WeightedEntry(constants::axes::d_type value, unsigned int count, double bin_center);

            /**
             * @brief Add the distance to this bin, and increase the counter by one.
             */
            void add(float distance, double value);

            /**
             * @brief Add twice the distance to this bin, and increase the counter by two.
             */
            void add2(float distance, double value);

            WeightedEntry operator+(const WeightedEntry& other) const;
            WeightedEntry& operator+=(const WeightedEntry& other);
            WeightedEntry operator-(const WeightedEntry& other) const;
            WeightedEntry& operator-=(const WeightedEntry& other);
            bool operator==(double other) const;
            
            constants::axes::d_type value = 0;
            unsigned int count = 0;
            double bin_center = 0;
        };

        WeightedEntry operator*(const WeightedEntry& entry, double factor);
        WeightedEntry operator*(double factor, const WeightedEntry& entry);
        std::ostream& operator<<(std::ostream& os, const WeightedEntry& entry);
    }
}