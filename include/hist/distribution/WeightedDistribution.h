#pragma once

#include <constants/Constants.h>

#include <vector>

namespace hist {
    namespace detail {
        struct Entry {
            int32_t count = 1; // Start at 1 to avoid division by zero.
            float content = 0;

            /**
             * @brief Add the distance to this bin, and increase the counter by one.
             */
            void add(float distance) {
                count++;
                content += distance;
            }
        };
    }
    static_assert(sizeof(detail::Entry) == 8, "hist::Entry is not 8 bytes long");

    struct WeightedDistribution {
        /**
         * @brief Reset the bin content tracker for the WeightedDistribution classes.
         */
        static void reset() {
            entries = std::vector<detail::Entry>(constants::axes::d_axis.bins);
        }

        /**
         * @brief Get a weighted distance axis based on the contents added to all WeightedDistribution instances since the last reset.
         */
        static std::vector<constants::axes::d_type> get_weighted_bins() {
            std::vector<constants::axes::d_type> weighted_bins(constants::axes::d_axis.bins, 0);
            std::transform(entries.begin(), entries.end(), constants::axes::d_vals.begin(), weighted_bins.begin(), [](const detail::Entry& entry, constants::axes::d_type d_val) {
                // if no counts were added to this bin, entry.count-1 = 0, which is negated to get 1
                // if any counts were added, entry.count-1 != 0, which is negated to get 0. thus the first term contributes only if the bin is empty
                return !(entry.count-1)*d_val + entry.content/(entry.count-1);
            });
            return weighted_bins;
        }

        inline static std::vector<detail::Entry> entries;
    };
}