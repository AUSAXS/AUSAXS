#pragma once

#include <container/Container3D.h>
#include <hist/distribution/WeightedEntry.h>

namespace hist {
    class WeightedDistribution3D {
        using type = constants::axes::d_type;

        public:
            WeightedDistribution3D(unsigned int size_x, unsigned int size_y, unsigned int size_z, type value);

            void add(unsigned int x, unsigned int y, float distance, type value);

            type& index(unsigned int x, unsigned int y, unsigned int z);
            const type& index(unsigned int x, unsigned int y, unsigned int z) const;

            const typename std::vector<detail::WeightedEntry>::const_iterator begin(unsigned int x, unsigned int y) const;
            const typename std::vector<detail::WeightedEntry>::const_iterator end(unsigned int x, unsigned int y) const;
            const typename std::vector<detail::WeightedEntry>::const_iterator begin() const;
            const typename std::vector<detail::WeightedEntry>::const_iterator end() const;

            typename std::vector<detail::WeightedEntry>::iterator begin(unsigned int x, unsigned int y);
            typename std::vector<detail::WeightedEntry>::iterator end(unsigned int x, unsigned int y);
            typename std::vector<detail::WeightedEntry>::iterator begin();
            typename std::vector<detail::WeightedEntry>::iterator end();

            std::size_t size_x() const;
            std::size_t size_y() const;
            std::size_t size_z() const;
            void resize(unsigned int size);

        private:
            container::Container3D<detail::WeightedEntry> data;
    };
}