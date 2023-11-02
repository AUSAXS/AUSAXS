#pragma once

#include <hist/distribution/WeightedEntry.h>
#include <container/Container1D.h>

namespace hist {
    class WeightedDistribution1D {
        using type = constants::axes::d_type;

        public:
            WeightedDistribution1D() = default;
            WeightedDistribution1D(unsigned int size, type value);

            void add(float distance, type value);

            type& index(unsigned int i);
            const type& index(unsigned int i) const;

            const typename std::vector<detail::WeightedEntry>::const_iterator begin() const;
            const typename std::vector<detail::WeightedEntry>::const_iterator end() const;

            typename std::vector<detail::WeightedEntry>::iterator begin();
            typename std::vector<detail::WeightedEntry>::iterator end();

            std::size_t size() const;
            bool empty() const;
            void resize(unsigned int size);

            /**
             * @brief Get a container with the values of the distribution.
             *        Note that this is a copy operation.
             * 
             * Complexity: O(n)
             */
            container::Container1D<type> get_container();

        private:
            container::Container1D<detail::WeightedEntry> data;
    };
}