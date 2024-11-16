#pragma once

#include <grid/GridFwd.h>
#include <math/MathFwd.h>
#include <utility/observer_ptr.h>

#include <vector>
#include <array>

namespace ausaxs::grid::detail {
    class RadialLineGenerator {
        public:
            RadialLineGenerator(observer_ptr<grid::Grid> grid, double radius, int divisions = 8);
            RadialLineGenerator(observer_ptr<grid::Grid> grid, std::array<double, 4> radii, int divisions = 8);
            virtual ~RadialLineGenerator();

            /**
             * @brief Generate the radial lines for the current grid. 
             *        This will only regenerate the radial lines if the grid width has changed. 
             */
            static void generate(std::array<double, 4> radius, int divisions);

        protected:
            static std::vector<Vector3<int>> rot_bins_1;
            static std::vector<Vector3<int>> rot_bins_2;
            static std::vector<Vector3<int>> rot_bins_3;
            static std::vector<Vector3<int>> rot_bins_4;
            static std::vector<Vector3<double>> rot_locs_abs;

        private:
            inline static double width = 0;
    };
}