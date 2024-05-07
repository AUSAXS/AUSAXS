#pragma once

#include <math/Vector3.h>
#include <crystal/CrystalFwd.h>
#include <grid/GridFwd.h>
#include <dataset/DatasetFwd.h>

#include <memory>
#include <vector>

namespace crystal {
    class CrystalScattering {
        public: 
            CrystalScattering(const grid::Grid& grid);

            CrystalScattering(const std::string& input);

            SimpleDataset calculate() const;
    
            SimpleDataset rotational_average(unsigned int n);

        private:
            std::shared_ptr<MillerGenerationStrategy> miller_strategy;

            /**
             * @brief Save a checkpoint file.
             */
            static void save_checkpoint(unsigned int start, const std::vector<Fval>& fvals);

            /**
             * @brief Load a checkpoint file.
             * 
             * @return The number of points loaded. Returns 0 if no checkpoint file exists.
             */
            static unsigned int load_checkpoint(std::vector<Fval>& fvals);

            void random_rotation();

            void rotate(const Vector3<double>& axis, double angle);

            void initialize();

            /**
             * @brief Convert the grid to individual points. 
             * 
             * All non-filled voxels are discarded. 
             */
            void convert_grid(const grid::Grid& grid) const; 
    };
}