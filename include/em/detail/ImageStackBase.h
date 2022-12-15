#pragma once

#include <vector>
#include <string>
#include <utility>

#include <em/Datatypes.h>
#include <em/Image.h>

namespace em {
    class PartialHistogramManager;

    /**
     * @brief \class ImageStackBase
     * 
     * A representation of a stack of images. 
     */
    class ImageStackBase {
        public:
            /**
             * @brief Constructor.
             * 
             * @param file Path to the input EM data file. 
             */
            ImageStackBase(std::string file);

            /**
             * @brief Constructor.
             * 
             * @param images The images for this stack.
             */
            ImageStackBase(const std::vector<Image>& images);

            /**
             * @brief Destructor.
             */
            ~ImageStackBase();

            /**
             * @brief Get a specific Image stored in this object. 
             * 
             * @param layer The vertical location of the Image. 
             */
            Image& image(unsigned int layer);

            /**
             * @brief Get a specific Image stored in this object. 
             * 
             * @param layer The vertical location of the Image. 
             */
            const Image& image(unsigned int layer) const;

            /**
             * @brief Prepare a ScatteringHistogram based on this object. 
             */
            hist::ScatteringHistogram get_histogram(double cutoff) const;

            /**
             * @brief Count the number of voxels for a given cutoff.
             */
            unsigned int count_voxels(double cutoff) const;

            /**
             * @brief Get the fitted ScatteringHistogram.
             */
            hist::ScatteringHistogram get_histogram(const std::shared_ptr<EMFit> res) const;

            /**
             * @brief Get the protein generated with the chosen cutoff value.
             */
            std::shared_ptr<Protein> get_protein(double cutoff) const;

            /**
             * @brief Get the header of the input file. 
             */
            std::shared_ptr<ccp4::Header> get_header() const;

            /**
             * @brief Get the number of images stored in this object.
             */
            unsigned int size() const;

            /**
             * @brief Get a reference to all images stored in this object. 
             */
            const std::vector<Image>& images() const;

            /**
             * @brief Save this structure as a .pdb file. 
             * 
             * @param cutoff The cutoff value. If positive, atoms will be generated at all pixel values higher than this. If negative, they will be generated at pixels lower than this. 
             * @param path Path to save location.
             */
            void save(double cutoff, std::string path) const;

            /**
             * @brief Get the limits on the q-values generated by discretizing the model scattering curve.
             */
            Limit get_limits() const;

            /**
             * @brief Get the mean density.
             */
            double mean() const;

            ObjectBounds3D minimum_volume(double cutoff);

            /**
             * @brief Get the cutoff corresponding to a PyMOL level. This is just the number of sigmas of the root-mean-square deviation.
             */
            double from_level(double sigma) const;

            /**
             * @brief Get the PyMOL level corresponding to a given cutoff. 
             */
            double to_level(double cutoff) const;

            /**
             * @brief Calculate the root-mean-square of this map. 
             */
            double rms() const;

            /**
             * @brief Get the histogram manager.
             */
            std::shared_ptr<em::PartialHistogramManager> get_histogram_manager() const;

            /**
             * @brief Determines the minimum bounds necessariy to describe the map for the given cutoff.
             * 
             * @param min_val The smallest possible value. Must be positive.
             * TODO Determine a better, more dynamic approach to determining this minimum. 
             */
            void set_minimum_bounds(double min_val);

        private:
            std::string filename;
            std::shared_ptr<ccp4::Header> header;               // The header of the map file. 
            std::vector<Image> data;                            // The actual image data. 
            unsigned int size_x, size_y, size_z;                // The number of pixels in each dimension.
            std::shared_ptr<em::PartialHistogramManager> phm;   // The histogram manager. Manages both the backing protein & its scattering curve. 
            
            void read(std::ifstream& istream, unsigned int byte_size);

            /**
             * @brief Get the data byte size of the CCP file. 
             */
            unsigned int get_byte_size() const;

            float& index(unsigned int x, unsigned int y, unsigned int z);

            float index(unsigned int x, unsigned int y, unsigned int z) const;
    };
}