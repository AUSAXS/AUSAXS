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

            std::vector<Vector3<int>> rot_bins_1;
            std::vector<Vector3<int>> rot_bins_2;
            std::vector<Vector3<int>> rot_bins_3;
            std::vector<Vector3<int>> rot_bins_4;
            std::vector<Vector3<double>> rot_locs_abs;

        private:
            /**
             * @brief Generate the radial lines for the current grid. 
             */
            void generate(double width, std::array<double, 4> radius, int divisions);
    };
}