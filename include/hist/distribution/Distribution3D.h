#pragma once

#include <container/Container3D.h>
#include <constants/Constants.h>

#include <cmath>

namespace hist {
    class Distribution3D {
        using type = constants::axes::d_type;

        public:
            Distribution3D() = default;
            Distribution3D(unsigned int size_x, unsigned int size_y, unsigned int size_z, type value);

            void add(unsigned int x, unsigned int y, float distance, type value);
            void add(unsigned int x, unsigned int y, int32_t i, type value);

            type& index(unsigned int x, unsigned int y, unsigned int z);
            const type& index(unsigned int x, unsigned int y, unsigned int z) const;

            const typename std::vector<type>::const_iterator begin(unsigned int x, unsigned int y) const;
            const typename std::vector<type>::const_iterator end(unsigned int x, unsigned int y) const;
            const typename std::vector<type>::const_iterator begin() const;
            const typename std::vector<type>::const_iterator end() const;

            typename std::vector<type>::iterator begin(unsigned int x, unsigned int y);
            typename std::vector<type>::iterator end(unsigned int x, unsigned int y);
            typename std::vector<type>::iterator begin();
            typename std::vector<type>::iterator end();

            std::size_t size_x() const;
            std::size_t size_y() const;
            std::size_t size_z() const;
            bool empty() const;
            void resize(unsigned int size);

            /**
             * @brief Get a container with the values of the distribution.
             *        Note that this is reference to the internal container.
             * 
             * Complexity: O(1)
             */
            container::Container3D<type>& get_container();

        private:
            container::Container3D<type> data;
    };
}