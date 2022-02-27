#pragma once

#include <fstream>
#include <vector>

#include <TH2D.h>

#include <em/datatypes.h>
#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <ScatteringHistogram.h>

using std::vector;

namespace em {
    class Image {
        public: 
            Image(std::shared_ptr<ccp4::Header> header);

            ~Image() = default;

            std::unique_ptr<TH2D> as_hist() const;

            void plot_atoms() const;

            void plot_without_solution() const;

            vector<Atom> generate_atoms(double cutoff) const;

            float index(unsigned int x, unsigned int y) const;
            float& index(unsigned int x, unsigned int y);

        private:
            std::shared_ptr<ccp4::Header> header;
            vector<vector<float>> data;
    };

    class ImageStack {
        public:
            ImageStack(string file);

            ImageStack(std::shared_ptr<ccp4::Header> header, std::ifstream& istream);

            ~ImageStack() = default;

            // 2D plot without background solution
            void plot_without_solution(unsigned int layer) const;

            // 3D plot
            void plot() const;

            void fit(string filename) const;

            Image& image(unsigned int layer);

            ScatteringHistogram calc_scattering_hist() const;

            std::unique_ptr<Grid> create_grid(double cutoff) const;

            std::unique_ptr<Protein> create_protein(double cutoff) const;

        private:
            std::shared_ptr<ccp4::Header> header;
            vector<Image> data;

            void read(std::ifstream& istream, size_t byte_size);

            size_t get_byte_size() const;

            float& index(unsigned int x, unsigned int y, unsigned int z);
            float index(unsigned int x, unsigned int y, unsigned int z) const;
    };
}