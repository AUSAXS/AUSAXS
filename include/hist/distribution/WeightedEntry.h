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

            /**
             * @brief This ensures that the WeightedEntry can be used with std::accumulate.
             */
            friend double operator+(const WeightedEntry& lhs, double v) {return lhs.value + v;}

            /**
             * @brief This ensures that the WeightedEntry can be used with std::accumulate.
             */
            friend double operator+(double v, const WeightedEntry& rhs) {return rhs.value + v;}

            friend WeightedEntry operator+(const WeightedEntry& lhs, const WeightedEntry& rhs) {
                WeightedEntry res = lhs;
                res += rhs;
                return res;
            }

            friend WeightedEntry operator*(const WeightedEntry& lhs, double rhs) {
                WeightedEntry res = lhs;
                res.value *= rhs;
                return res;
            }

            type distance;
            type value;
            type count;
    };
}