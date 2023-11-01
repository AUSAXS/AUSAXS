#pragma once

#include <container/Container2D.h>
#include <constants/Constants.h>

namespace hist {
    class Distribution2D {
        using type = constants::axes::d_type;
        public:
            Distribution2D() = default;
            Distribution2D(unsigned int size_x, unsigned int size_y, type value);

            /**
             * @brief Add a value for a given distance.
             */
            void add(unsigned int x, float distance, type value);
            void add(unsigned int x, int32_t i, type value);

            type& index(unsigned int x, unsigned int y);
            const type& index(unsigned int x, unsigned int y) const;

            const typename std::vector<type>::const_iterator begin(unsigned int x) const;
            const typename std::vector<type>::const_iterator end(unsigned int x) const;
            const typename std::vector<type>::const_iterator begin() const;
            const typename std::vector<type>::const_iterator end() const;

            typename std::vector<type>::iterator begin(unsigned int x);
            typename std::vector<type>::iterator end(unsigned int x);
            typename std::vector<type>::iterator begin();
            typename std::vector<type>::iterator end();

            std::size_t size_x() const;
            std::size_t size_y() const;
            void resize(unsigned int size);

        private:
            container::Container2D<type> data;
    };
}