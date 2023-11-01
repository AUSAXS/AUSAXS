#pragma once

#include <container/Container2D.h>
#include <hist/distribution/WeightedEntry.h>

namespace hist {
    template<typename T>
    class WeightedDistribution2D {
        public:
            WeightedDistribution2D(unsigned int size_x, unsigned int size_y, T value) : data(size_x, size_y, value) {}

            void add(unsigned int x, float distance, T value) {
                data.index(x, std::round(distance)).add(distance, value);
            }

            const T index(unsigned int x, unsigned int y) const {return data.index(x, y).value;}

            const typename std::vector<T>::const_iterator begin(unsigned int x) const {return data.begin(x);}
            const typename std::vector<T>::const_iterator end(unsigned int x) const {return data.end(x);}
            const typename std::vector<T>::const_iterator begin() const {return data.begin();}
            const typename std::vector<T>::const_iterator end() const {return data.end();}

            typename std::vector<T>::iterator begin(unsigned int x) {return data.begin(x);}
            typename std::vector<T>::iterator end(unsigned int x) {return data.end(x);}
            typename std::vector<T>::iterator begin() {return data.begin();}
            typename std::vector<T>::iterator end() {return data.end();}

            std::size_t size_x() const {return data.size_x();}
            std::size_t size_y() const {return data.size_y();}

        private:
            container::Container2D<detail::WeightedEntry<T>> data;
    };
}