#pragma once

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

            /**
             * @brief Get the mean density. 
             */
            double mean() const;

            /**
             * @brief Get the minimum and maximum density.
             */
            Limit limits() const;

            float index(unsigned int x, unsigned int y) const;
            float& index(unsigned int x, unsigned int y);

        private:
            std::shared_ptr<ccp4::Header> header;
            vector<vector<float>> data;
            unsigned int z;
    };
}