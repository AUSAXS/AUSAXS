#pragma once

#include <constants/Constants.h>
#include <container/ThreadLocalWrapper.h>

#include <vector>

namespace hist {
    namespace detail {
        struct Entry {
            int count = 0;
            double content = 0;

            /**
             * @brief Add the distance to this bin, and increase the counter by one.
             */
            void add(double distance) {
                ++count;
                content += distance;
            }
        };

        // struct ThreadContainer;
        // inline static std::vector<ThreadContainer*> thread_containers;
        // struct ThreadContainer : std::vector<Entry> {
        //     ThreadContainer() {
        //         entries = std::vector<Entry>(constants::axes::d_axis.bins);
        //         thread_containers.push_back(this);
        //         std::cout << "thread_containers.size() = " << thread_containers.size() << std::endl;
        //     }

        //     Entry& operator[](unsigned int i) {
        //         return entries[i];
        //     }

        //     const Entry& operator[](unsigned int i) const {
        //         return entries[i];
        //     }

        //     std::size_t size() const {
        //         return entries.size();
        //     }

        //     std::vector<Entry> entries;
        // };
    }

    struct WeightedDistribution {
        /**
         * @brief Reset the bin content tracker for the WeightedDistribution classes.
         */
        static void reset() {
            // for (unsigned int i = 0; i < detail::thread_containers.size(); ++i) {
            //     detail::thread_containers[i]->entries = std::vector<detail::Entry>(constants::axes::d_axis.bins);
            // }
            entries = std::vector<detail::Entry>(constants::axes::d_axis.bins);
        }

        /**
         * @brief Get a weighted distance axis based on the contents added to all WeightedDistribution instances since the last reset.
         */
        static std::vector<constants::axes::d_type> get_weighted_bins() {
            std::vector<constants::axes::d_type> weighted_bins(constants::axes::d_axis.bins, 0);

            // merge the thread_local containers
            std::vector<detail::Entry> merged(constants::axes::d_axis.bins);
            for (auto& entry : entries.get_all()) {
                for (unsigned int i = 0; i < constants::axes::d_axis.bins; ++i) {
                    merged[i].count += entry.get()[i].count;
                    merged[i].content += entry.get()[i].content;
                }
            }

            // for (unsigned int i = 0; i < detail::thread_containers.size(); ++i) {
            //     for (unsigned int j = 0; j < constants::axes::d_axis.bins; ++j) {
            //         entries[j].count += detail::thread_containers[i]->entries[j].count;
            //         entries[j].content += detail::thread_containers[i]->entries[j].content;
            //     }
            // }

            // calculate the weighted bins
            std::transform(merged.begin(), merged.end(), constants::axes::d_vals.begin(), weighted_bins.begin(), [](const detail::Entry& entry, constants::axes::d_type d_val) -> double {
                // if no counts were added to this bin, entry.count == 0, which when negated is 1
                // if any counts were added, entry.count != 0, which when negated is 0. 
                // thus the first term is non-zero iff no counts were added to this bin, and similarly a 1 is added to the divisor of the second term iff no counts were added
                // we do this to avoid a nested if statement
                return (!entry.count)*d_val + entry.content/(entry.count + !entry.count);
            });
            weighted_bins[0] = 0; // the first bin must always be 0
            return weighted_bins;
        }

        inline static container::ThreadLocalWrapper<std::vector<detail::Entry>> entries = std::vector<detail::Entry>(constants::axes::d_axis.bins);
        // inline static thread_local detail::ThreadContainer entries;
    };
}