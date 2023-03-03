#pragma once

#include <hydrate/Grid.h>
#include <crystal/miller/MillerGenerationStrategy.h>
#include <dataset/SimpleDataset.h>

#include <memory>

namespace crystal {
    class CrystalScattering {
        public: 
            CrystalScattering(const Grid& grid);

            CrystalScattering(const std::string& input);

            SimpleDataset calculate() const;

        private:
            std::shared_ptr<MillerGenerationStrategy> miller_strategy;

            void initialize();

            /**
             * @brief Convert the grid to individual points. 
             * 
             * All non-filled voxels are discarded. 
             */
            void convert_grid(const Grid& grid) const; 
    };
}