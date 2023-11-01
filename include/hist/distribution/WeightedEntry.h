#pragma once

#include <constants/Constants.h>

namespace hist::detail {
    /**
     * @brief A small struct which keeps track of which distances has been added to it. 
     *        typehis is useful for later using weighted bins.
     */
    class WeightedEntry {
        using type = constants::axes::d_type;

        public:
            WeightedEntry() = default;
            WeightedEntry(type value) : distance(0), value(value) {}

            void add(float distance, type value) {
                distance += distance;
                value += value;
                count++;
            }

            WeightedEntry& operator+=(const WeightedEntry& rhs) {
                distance += rhs.distance;
                value += rhs.value;
                count += rhs.count;
                return *this;
            }

            WeightedEntry& operator-=(const WeightedEntry& rhs) {
                distance -= rhs.distance;
                value -= rhs.value;
                count -= rhs.count;
                return *this;
            }

            type distance;
            type value;
            type count;
    };
}