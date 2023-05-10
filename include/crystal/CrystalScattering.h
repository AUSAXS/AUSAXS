#pragma once

#include <memory>

namespace grid {class Grid;}
class SimpleDataset;
namespace crystal {
    class MillerGenerationStrategy;
    class CrystalScattering {
        public: 
            CrystalScattering(const grid::Grid& grid);

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
            void convert_grid(const grid::Grid& grid) const; 
    };
}