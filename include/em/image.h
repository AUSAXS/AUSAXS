#pragma once

#include <fstream>
#include <vector>
#include <list>

#include <TH2D.h>

#include <em/datatypes.h>
#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <ScatteringHistogram.h>

using std::vector, std::list;

namespace em {
    /**
     * @brief \class Image
     * 
     * Supporting class for ImageStack. This is not meant to be instantiated elsewhere. 
     */
    class Image {
        public: 
            /**
             * @brief Constructor.
             * 
             * @param header Header of the parent ImageStack. 
             * @param layer The layer number of this Image. 
             */
            Image(std::shared_ptr<ccp4::Header> header, unsigned int layer = 0);

            ~Image() = default;

            std::unique_ptr<TH2D> as_hist() const;

            list<Atom> generate_atoms(double cutoff) const;

            /**
             * @brief Set the z location of this object. 
             */
            void set_z(unsigned int z);

            float index(unsigned int x, unsigned int y) const;
            float& index(unsigned int x, unsigned int y);

        private:
            std::shared_ptr<ccp4::Header> header;
            vector<vector<float>> data;
            unsigned int z;
    };

    /**
     * @brief \class ImageStack
     * 
     * A representation of a stack of images. 
     */
    class ImageStack {
        public:
            /**
             * @brief Constructor.
             * 
             * @param file Path to the input EM data file. 
             */
            ImageStack(string file);

            /**
             * @brief Destructor.
             */
            ~ImageStack() = default;

            /**
             * @brief Fit the cutoff value with the input experimental data file. 
             * 
             * @param filename Path to the measurement file. 
             */
            void fit(string filename) const;

            /**
             * @brief Get a specific Image stored in this object. 
             * 
             * @param layer The vertical location of the Image. 
             */
            Image& image(unsigned int layer);

            /**
             * @brief Prepare a ScatteringHistogram based on this object. 
             */
            ScatteringHistogram calc_scattering_hist() const;

            /**
             * @brief Create a new Grid based on this object. 
             * 
             * @param cutoff The cutoff value. If positive, atoms will be generated at all pixel values higher than this. If negative, they will be generated at pixels lower than this. 
             */
            std::unique_ptr<Grid> create_grid(double cutoff) const;

            /**
             * @brief Create a new Protein based on this object. 
             * 
             * @param cutoff The cutoff value. If positive, atoms will be generated at all pixel values higher than this. If negative, they will be generated at pixels lower than this. 
             */
            std::unique_ptr<Protein> create_protein(double cutoff) const;

            /**
             * @brief Get the header of the input file. 
             */
            std::shared_ptr<ccp4::Header> get_header() const;

            /**
             * @brief Get the number of images stored in this object.
             */
            size_t size() const;

            /**
             * @brief Get a reference to all images stored in this object. 
             */
            const vector<Image>& images() const;

            /**
             * @brief Save this structure as a .pdb file. 
             * 
             * @param path Path to save location.
             * @param cutoff The cutoff value. If positive, atoms will be generated at all pixel values higher than this. If negative, they will be generated at pixels lower than this. 
             */
            void save(string path, double cutoff) const;

        private:
            std::shared_ptr<ccp4::Header> header;
            vector<Image> data;

            void read(std::ifstream& istream, size_t byte_size);

            size_t get_byte_size() const;

            float& index(unsigned int x, unsigned int y, unsigned int z);
            float index(unsigned int x, unsigned int y, unsigned int z) const;
    };
}