#pragma once

#include <container/Container3D.h>
#include <hist/distribution/WeightedEntry.h>

namespace hist {
    template<typename T>
    class WeightedDistribution3D {
        public:
            WeightedDistribution3D(unsigned int size_x, unsigned int size_y, unsigned int size_z, T value) : data(size_x, size_y, size_z, value) {}

            void add(unsigned int x, float distance, T value) {
                data.index(x, std::round(distance)).add(distance, value);
            }

            const T index(unsigned int x, unsigned int y) const {return data.index(x, y).value;}

            const typename std::vector<T>::const_iterator begin(unsigned int x, unsigned int y) const {return data.begin(x, y);}
            const typename std::vector<T>::const_iterator end(unsigned int x, unsigned int y) const {return data.end(x, y);}
            const typename std::vector<T>::const_iterator begin() const {return data.begin();}
            const typename std::vector<T>::const_iterator end() const {return data.end();}

            typename std::vector<T>::iterator begin(unsigned int x, unsigned int y) {return data.begin(x, y);}
            typename std::vector<T>::iterator end(unsigned int x, unsigned int y) {return data.end(x, y);}
            typename std::vector<T>::iterator begin() {return data.begin();}
            typename std::vector<T>::iterator end() {return data.end();}

            std::size_t size_x() const {return data.size_x();}
            std::size_t size_y() const {return data.size_y();}
            std::size_t size_z() const {return data.size_z();}

        private:
            container::Container3D<detail::WeightedEntry<T>> data;
    };
}