// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <em/ObjectBounds2D.h>
#include <data/DataFwd.h>
#include <em/detail/header/HeaderFwd.h>
#include <em/detail/EMAtom.h>
#include <hist/HistFwd.h>
#include <utility/observer_ptr.h>

#include <list>

namespace ausaxs::em {
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
            Image(observer_ptr<em::detail::header::IMapHeader> header, unsigned int layer = 0);

            Image(const Matrix<float>& data);

            Image(const Matrix<float>& data, observer_ptr<em::detail::header::IMapHeader> header, unsigned int layer);
            
            ~Image() = default;

            hist::Histogram2D as_hist() const;

            std::list<data::EMAtom> generate_atoms(double cutoff) const;

            /**
             * @brief Count the number of voxels larger than a given cutoff value.
             */
            unsigned int count_voxels(double cutoff) const;

            /**
             * @brief Set the z location of this object. 
             */
            void set_z(unsigned int z);

            /**
             * @brief Get the z location of this object. 
             */
            unsigned int get_z() const;

            /**
             * @brief Get the mean density. 
             */
            double mean() const;

            /**
             * @brief Get the minimum and maximum density.
             */
            Limit limits() const;

            /**
             * @brief Get the current bounds of this image. 
             */
            const ObjectBounds2D& get_bounds() const;

            /**
             * @brief Set the header. 
             */
            void set_header(observer_ptr<em::detail::header::IMapHeader> header);

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

            /**
             * @brief Get the squared sum of all entries in this image. 
             *        This is intended as a helper method for the RMS method of the ImageStack class.
             */
            double squared_sum() const;

            bool operator==(const Image& other) const;

            /**
             * @brief Get a string representation of this object. 
             */
            std::string to_string() const;

            /**
             * @brief Get the data of this image. 
             */
            const Matrix<float>& get_data() const;

            unsigned int N; // The number of rows.  
            unsigned int M; // The number of columns.
        private:
            observer_ptr<em::detail::header::IMapHeader> header; // Pointer to the header of the parent ImageStack.
            Matrix<float> data;                                  // The actual data storage. 
            unsigned int z;                                      // The z-index of this image in the ImageStack. 
            ObjectBounds2D bounds;
    };
}