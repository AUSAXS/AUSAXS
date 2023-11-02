#pragma once

#include <container/Container2D.h>
#include <hist/distribution/WeightedEntry.h>

namespace hist {
    class WeightedDistribution2D {
        using type = constants::axes::d_type;

        public:
            WeightedDistribution2D() = default;
            WeightedDistribution2D(unsigned int size_x, unsigned int size_y, type value);

            void add(unsigned int x, float distance, type value);

            type& index(unsigned int x, unsigned int y);
            const type& index(unsigned int x, unsigned int y) const;

            const typename std::vector<detail::WeightedEntry>::const_iterator begin(unsigned int x) const;
            const typename std::vector<detail::WeightedEntry>::const_iterator end(unsigned int x) const;
            const typename std::vector<detail::WeightedEntry>::const_iterator begin() const;
            const typename std::vector<detail::WeightedEntry>::const_iterator end() const;

            typename std::vector<detail::WeightedEntry>::iterator begin(unsigned int x);
            typename std::vector<detail::WeightedEntry>::iterator end(unsigned int x);
            typename std::vector<detail::WeightedEntry>::iterator begin();
            typename std::vector<detail::WeightedEntry>::iterator end();

            std::size_t size_x() const;
            std::size_t size_y() const;
            bool empty() const;
            void resize(unsigned int size);

            /**
             * @brief Get a container with the values of the distribution.
             *        Note that this is a copy operation.
             * 
             * Complexity: O(n)
             */
            container::Container2D<type> get_container();

        private:
            container::Container2D<detail::WeightedEntry> data;
    };
}