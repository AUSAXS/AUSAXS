#pragma once

#include <grid/GridFwd.h>
#include <math/MathFwd.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace grid::detail {
    class RadialLineGenerator {
        public:
            RadialLineGenerator(observer_ptr<grid::Grid> grid);
            ~RadialLineGenerator() = default;

            /**
             * @brief Generate the radial lines for the current grid. 
             *        This will only regenerate the radial lines if the grid width has changed. 
             */
            static void generate(int divisions = 8);

        protected:
            static std::vector<Vector3<int>> rot_bins_1rh; // rotation bins at 1rh radius
            static std::vector<Vector3<int>> rot_bins_3rh; // rotation bins at 3rh radius
            static std::vector<Vector3<int>> rot_bins_5rh; // rotation bins at 5rh radius
            static std::vector<Vector3<int>> rot_bins_7rh; // rotation bins at 7rh radius
            static std::vector<Vector3<double>> rot_locs;  // absolute locations of the rotation bins

        private:
            inline static double width = 0;
    };
}