#pragma once

#include <fstream>
#include <vector>

#include <em/datatypes.h>
#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <ScatteringHistogram.h>

using std::vector;

namespace em {
    class ImageStack {
        public:
            ImageStack(string file);

            ImageStack(const ccp4::Header& header, std::ifstream& istream);

            ~ImageStack() = default;

            void plot_no_solution(unsigned int layer) const;

            void plot(unsigned int layer) const;

            void plot3d() const;

            void fit(string filename) const;

        private:
            ccp4::Header header;
            // vector<T> data;
            vector<vector<vector<float>>> data;

            void read(std::ifstream& istream, size_t byte_size);

            ScatteringHistogram calc_scattering_hist() const;

            std::unique_ptr<Grid> create_grid(double cutoff) const;

            std::unique_ptr<Protein> create_protein(double cutoff) const;

            size_t get_byte_size() const;

            float& index(unsigned int x, unsigned int y, unsigned int z);
            float index(unsigned int x, unsigned int y, unsigned int z) const;
    };
}