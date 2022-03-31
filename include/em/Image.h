#pragma once

#include <vector>
#include <list>

#include <TH2D.h>

#include <em/datatypes.h>
#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <ScatteringHistogram.h>
#include <math/Matrix.h>

/**
 * @brief Describes the bounds of some object contained within a 2D matrix. 
 */
class MatrixBounds {
    public:
        MatrixBounds(unsigned int size) : bounds(size) {}

        Limit& get_bounds(unsigned int x) {return bounds[x];}
        
        void set_bounds(unsigned int x, const Limit& bounds) {this->bounds[x] = bounds;}

        Limit& operator[](unsigned int x) {return bounds[x];}

        const Limit& operator[](unsigned int x) const {return bounds[x];}

        unsigned int size() const {return bounds.size();}

        vector<Limit> bounds;
};

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

            Image(const Matrix<float>& data);

            ~Image() = default;

            std::unique_ptr<TH2D> as_hist() const;

            std::list<Atom> generate_atoms(double cutoff) const;

            /**
             * @brief Set the z location of this object. 
             */
            void set_z(unsigned int z);

            /**
             * @brief Get the mean density. 
             */
            double mean() const;

            /**
             * @brief Get the minimum and maximum density.
             */
            Limit limits() const;

            /**
             * @brief Get the minimum area covering all pixels with a density more extreme than the cutoff. 
             * 
             * @param cutoff The density cutoff value. If positive, the area will cover pixels with a density @a higher than this. 
             *               If negative, it will cover pixels with density @a lower than this. 
             * 
             * @return A vector containing the minimum and maximum x-indices covering the area. 
             */
            MatrixBounds minimum_area(double cutoff) const;

            float index(unsigned int x, unsigned int y) const;
            float& index(unsigned int x, unsigned int y);

        private:
            std::shared_ptr<ccp4::Header> header;
            unsigned int N; // The number of rows.  
            unsigned int M; // The number of columns.
            Matrix<float> data; // The actual data storage. 
            unsigned int z; // The z-index of this image in the ImageStack. 
    };
}