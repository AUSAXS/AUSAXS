#pragma once

#include <container/Container1D.h>
#include <constants/Constants.h>

#include <cmath>

namespace hist {
    class Distribution1D {
        using type = constants::axes::d_type;

        public:
            Distribution1D() = default;
            Distribution1D(unsigned int size, type value);

            /**
             * @brief Add a value for a given distance.
             */
            void add(float distance, type value);
            void add(int32_t i, type value);

            type& index(unsigned int i);
            const type& index(unsigned int i) const;

            const typename std::vector<type>::const_iterator begin() const;
            const typename std::vector<type>::const_iterator end() const;

            typename std::vector<type>::iterator begin();
            typename std::vector<type>::iterator end();

            std::size_t size() const;
            void resize(unsigned int size);

            std::vector<type>& get_counts();

        private:
            container::Container1D<type> data;
    };
}