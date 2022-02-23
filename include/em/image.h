#pragma once

#include <fstream>
#include <vector>

#include <em/datatypes.h>

using std::vector;

namespace em {
    class Image {
        public:
            Image(const ccp4::Header& header, std::ifstream& istream);

            void plot_no_solution(unsigned int layer) const;

            void plot(unsigned int layer) const;

            void plot3d() const;

        private:
            const ccp4::Header header;
            // vector<T> data;
            vector<vector<vector<float>>> data;

            void read(std::ifstream& istream, size_t byte_size);

            size_t get_byte_size() const;

            float index(unsigned int x, unsigned int y, unsigned int z) const;
    };
}