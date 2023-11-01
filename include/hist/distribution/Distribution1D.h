#pragma once

#include <container/Container1D.h>

namespace hist {
    template<typename T>
    class Distribution1D {
        public:
            Distribution1D(unsigned int size, T value) : data(size, value) {}

            /**
             * @brief Add a value for a given distance.
             */
            void add(float distance, T value) {data.index(std::round(distance)) += value;}
            void add(unsigned int i, T value) {data.index(i) += value;}

            const T index(unsigned int i) const {return data.index(i);}

            const typename std::vector<T>::const_iterator begin() const {return data.begin();}
            const typename std::vector<T>::const_iterator end() const {return data.end();}

            typename std::vector<T>::iterator begin() {return data.begin();}
            typename std::vector<T>::iterator end() {return data.end();}

            std::size_t size() const {return data.size();}

        private:
            container::Container1D<T> data;
    };
}