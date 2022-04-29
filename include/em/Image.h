#pragma once

#include <vector>
#include <list>

#include <TH2D.h>

#include <em/Datatypes.h>
#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <ScatteringHistogram.h>
#include <math/Matrix.h>

/**
 * @brief Describes the bounds of some object contained within a 2D matrix. 
 */
class ObjectBounds2D {
    public:
        ObjectBounds2D(unsigned int size_x, unsigned int size_y) : bounds(size_x, Limit(0, size_y)), N(size_x), M(size_y) {}

        Limit& operator[](unsigned int x) {return bounds[x];}

        const Limit& operator[](unsigned int x) const {return bounds[x];}

        unsigned int size() const {return bounds.size();}

        bool empty() const {return bounds.empty();}

        unsigned int bounded_area() const {return std::accumulate(bounds.begin(), bounds.end(), 0, [] (unsigned int area, const Limit& limit) {return area += limit.max+1 - limit.min;});}

        unsigned int total_area() const {return N*M;}

    private:
        std::vector<Limit> bounds;
        unsigned int N, M;
};

class ObjectBounds3D {
    public: 
        ObjectBounds3D(unsigned int size_x, unsigned int size_y, unsigned int size_z) : bounds(size_z, ObjectBounds2D(size_x, size_y)), size_x(size_x), size_y(size_y), size_z(size_z) {}

        ObjectBounds2D& operator[](unsigned int z) {return bounds[z];}

        const ObjectBounds2D& operator[](unsigned int z) const {return bounds[z];}

        unsigned int total_volume() const {return size_x*size_y*size_z;}

        unsigned int bounded_volume() const {return std::accumulate(bounds.begin(), bounds.end(), 0, [] (unsigned int volume, const ObjectBounds2D& bound) {return volume += bound.bounded_area();});}

    private:
        std::vector<ObjectBounds2D> bounds;
        unsigned int size_x, size_y, size_z;
};

namespace em {
    /**
     * @brief Supporting class for ImageStack. This is not meant to be instantiated elsewhere. 
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

            Image(const Matrix<float>& data, std::shared_ptr<ccp4::Header> header, unsigned int layer);
            
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

            const ObjectBounds2D& get_bounds() const;

            /**
             * @brief Get the minimum area covering all pixels with a density more extreme than the cutoff. 
             * 
             * @param cutoff The density cutoff value. If positive, the area will cover pixels with a density @a higher than this. 
             *               If negative, it will cover pixels with density @a lower than this. 
             * 
             * @return A vector containing the minimum and maximum x-indices covering the area. 
             */
            const ObjectBounds2D& setup_bounds(double cutoff);

            float index(unsigned int x, unsigned int y) const;
            float& index(unsigned int x, unsigned int y);

            unsigned int N; // The number of rows.  
            unsigned int M; // The number of columns.
        private:
            std::shared_ptr<ccp4::Header> header;
            Matrix<float> data; // The actual data storage. 
            unsigned int z; // The z-index of this image in the ImageStack. 
            ObjectBounds2D bounds;
    };
}