#pragma once

#include <constants/Constants.h>

#include <vector>
#include <atomic>

#include <iostream> //! remove
namespace hist {
    namespace detail {
        struct Entry {
            std::atomic<int> count = 0;
            std::atomic<float> content = 0;

            /**
             * @brief Add the distance to this bin, and increase the counter by one.
             */
            void add(float distance) {
                ++count;
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
                // if no counts were added to this bin, entry.count == 0, which when negated is 1
                // if any counts were added, entry.count != 0, which when negated is 0. 
                // thus the first term is non-zero iff no counts were added to this bin, and similarly a 1 is added to the divisor of the second term iff no counts were added
                // we do this to avoid a nested if statement
                return (!entry.count)*d_val + entry.content/(entry.count + !entry.count);
            });

            for (unsigned int i = 0; i < 20; ++i) {
                std::cout << "entry[" << i << "]: " << entries[i].content << " / " << entries[i].count << " = " << entries[i].content/entries[i].count << " = " << weighted_bins[i] << std::endl;
            }

            return weighted_bins;
        }

        inline static std::vector<detail::Entry> entries;
    };
}